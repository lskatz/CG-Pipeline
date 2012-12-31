#!/usr/bin/env perl

# Author: Lee Katz <lkatz@cdc.gov>
# create a genbank file from the annotation sql files

use warnings;
use strict;
use FindBin;
use lib "$FindBin::RealBin/../lib";
use AKUtils qw/logmsg/;
use Bio::Seq;
use Bio::Seq::RichSeq;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Bio::Annotation::Reference;
use File::Basename;
use Data::Dumper;
use Getopt::Long;

exit(&main());

sub main(){
  my $settings={appname=>'cgpipeline'};
  $settings=AKUtils::loadConfig($settings);
  die usage() if(@ARGV < 1);
  GetOptions($settings,qw(organism=s prediction=s inputdir=s gb=s gff=s)) or die;
  $$settings{numcpus}||=AKUtils::getNumCPUs();
  my $annotationDir=$$settings{inputdir};
  die "ERROR: cannot read annotation directory at $annotationDir" if(!-d $annotationDir || !-x $annotationDir);
  my $predictionGbk=$$settings{prediction};
  die "ERROR: cannot read prediction file at $predictionGbk" if(!-f $predictionGbk);
  my $outfile=$$settings{gb};
  die "ERROR: need genbank parameter".usage() if(!$outfile);
  
  logmsg "Putting all annotations in $annotationDir into memory";
  my $annotationThread=threads->new(\&readAnnotationDir,$annotationDir,$settings);
  
  logmsg "Transforming the prediction genbank into a set of features";
  my $gbkToFeaturesThread=threads->new(\&gbkToFeatures,$predictionGbk,$settings);
  
  # join threads before continuing
  my $tmp=$gbkToFeaturesThread->join;
  my($predfeatures,$metaData)=@$tmp;
  my $annotation=$annotationThread->join;
  
  logmsg "Combining the features from the prediction file and the annotations";
  combinePredictionAndAnnotation($outfile,$predfeatures,$metaData,$annotation,$settings);
  
  return 0;
}

sub combinePredictionAndAnnotation{
  my ($outfile,$predfeatures,$metaData,$annotation,$settings)=@_;
  my %annotationDirMap=annotationFieldMap();
  my @annotationKeys=keys(%annotationDirMap);
  
  # use uniprot and cogs information
  my %uniprotEvidence=uniprotEvidence($$settings{inputdir},\%annotationDirMap,$settings);
  my %uniprot=uniprotSql($$settings{inputdir},\%annotationDirMap,$settings);
  my %prot2cogs=prot2cogMapping($settings);
  my $extraAnnotationInfo={uniprotEvidence=>\%uniprotEvidence,uniprot=>\%uniprot,prot2cogs=>\%prot2cogs};
  
  logmsg "One contig/chromosome per dot: ";
  my $out=Bio::SeqIO->new(-file=>">$outfile",-format=>"genbank");
  while(my ($seqid,$meta)=each(%$metaData)){
    $|++;print ".";$|--; # print out a dot per contig and flush the buffer
    my $sourceFeat=$$meta{source};
    my $seq=$sourceFeat->seq;
    $seq=Bio::Seq->new(-seq=>$seq->seq,-id=>$seq->id);
    $seq->add_SeqFeature($sourceFeat);
    $seq->annotation->add_Annotation('reference',$$metaData{$seqid}{'reference'});
    my @sortedLocus=sort{
      return $$predfeatures{$seqid}{'gene'}{$a} <=> $$predfeatures{$seqid}{'gene'}{$b};
    } keys(%{$$predfeatures{$seqid}{'gene'}});
    
    # add annotations from the annotations directory
    # TODO this can be multithreaded
    for my $locus_tag(@sortedLocus){
      my $locusAnnotation=$$annotation{$locus_tag};
      for(@annotationKeys){
        $$locusAnnotation{$_}=[] if(!$$locusAnnotation{$_})
      }
      
      my @geneFeat; # the following subroutines should internally sort this correctly
      my $geneFeat=$$predfeatures{$seqid}{'gene'}{$locus_tag}->clone;
      if($$predfeatures{$seqid}{'CDS'}{$locus_tag}){
        @geneFeat=interpretCdsFeat($geneFeat,$$predfeatures{$seqid}{'CDS'}{$locus_tag},$locusAnnotation,$extraAnnotationInfo,$settings);
      } elsif($$predfeatures{$seqid}{'rRNA'}{$locus_tag}){
        @geneFeat=interpretRRnaFeat($geneFeat,$$predfeatures{$seqid}{'rRNA'}{$locus_tag},$settings);
      } elsif($$predfeatures{$seqid}{'tRNA'}{$locus_tag}){
        @geneFeat=interpretTRnaFeat($geneFeat,$$predfeatures{$seqid}{'tRNA'}{$locus_tag},$settings);
      } else {
        logmsg "WARNING: I could not understand what kind of gene $locus_tag is. Skipping annotation.";
        @geneFeat=($geneFeat);
      }
      
      $seq->add_SeqFeature(@geneFeat);
    }
    
    $out->write_seq($seq);
  }
  $out->close;
  print "\n"; # newline after all the "update dots"
  
  return 1;
}

sub interpretCdsFeat{
  my($geneFeat,$cdsFeatOrig,$locusAnnotation,$extraAnnotationInfo,$settings)=@_;
  my @feat;
  my(@miscFeat,@otherFeat);
  my $cdsFeat=$cdsFeatOrig->clone;
  # modify the gene and cds features by reference
  interpretUniprot($geneFeat,$cdsFeat,$locusAnnotation,$settings);
  
  # other whole-gene annotators: IS elements, VFs, ...
  my $isFeat=interpretIs($cdsFeat,$locusAnnotation,$settings);
  push(@miscFeat,$isFeat) if($isFeat);
  my $vfFeat=interpretVf($cdsFeat,$locusAnnotation,$settings);
  push(@miscFeat,$vfFeat) if($vfFeat);
  my $cogsFeat=interpretCogs($cdsFeat,$locusAnnotation,$extraAnnotationInfo,$settings);
  push(@miscFeat,$cogsFeat) if($cogsFeat);
  # TODO PDB
  
  ## annotations on portions of the gene
  
  # signal peptide
  my $signalpFeat=interpretSignalP($cdsFeat,$locusAnnotation,$settings);
  push(@miscFeat,$signalpFeat) if($signalpFeat); # $signalpFeat is 0 if not present
      
  # transmembrane helix
  my $tmFeat=interpretTmhmm($cdsFeat,$locusAnnotation,$settings);
  push(@miscFeat,$tmFeat) if($tmFeat);
      
  # domain hits
  my @iprFeat=interpretIpr($cdsFeat,$$locusAnnotation{ipr_matches},$settings);
  push(@otherFeat,@iprFeat);
  @otherFeat=sort {$a->start<=>$b->start} @otherFeat;
  
  push(@feat,$geneFeat,$cdsFeat,@miscFeat,@otherFeat);
  return @feat;
}

# I think that rRNAs are already annotated well
sub interpretRRnaFeat{
  my($geneFeat,$rnaFeatOrig,$locusAnnotation,$settings)=@_;
  return($geneFeat,$rnaFeatOrig);
}
# I think that tRNAs are already annotated well
sub interpretTRnaFeat{
  my($geneFeat,$rnaFeatOrig,$locusAnnotation,$settings)=@_;
  return($geneFeat,$rnaFeatOrig);
}

sub gbkToFeatures{
  my($gbk,$settings)=@_;
  my %feat;
  my %metaData;
  # reference information
  my $reference="";
  if($$settings{'pipeline_reference'}){
    my ($gb_authors,$gb_title,$gb_journal,$gb_pubmed)=readSqlLine($$settings{'pipeline_reference'},$settings);
    $reference=Bio::Annotation::Reference->new(-title=>$gb_title,-authors=>$gb_authors,-location=>$gb_journal,-pubmed=>$gb_pubmed) if(defined($$settings{'pipeline_reference'}));
  } else {
    logmsg "WARNING: could not find the parameter pipeline_reference in your config file. The genbank reference line will not be filled in."
  }
    
  my $in=Bio::SeqIO->new(-file=>$gbk);
  while(my $seq=$in->next_seq){
    my @feats = $seq->all_SeqFeatures();
    
    # looking at gene features only
    for my $feat(@feats){
      #delete($$feat{'_gsf_seq'}); # save memory and simplify things by removing the embedded seq object (hopefully not buggy)
      #next if($feat->primary_tag()!~/gene/i); # only look at genes
      my $primary_tag=$feat->primary_tag();
      next if(!$feat->has_tag('locus_tag'));
      my @locus_tag=$feat->get_tag_values('locus_tag');
      die "ERROR: a locus in $gbk has more than one locus_tag" if(@locus_tag>1);
      my $locus_tag=$locus_tag[0];
      die "ERROR: locus tag $locus_tag exists twice." if($feat{$seq->id}{$primary_tag}{$locus_tag});
      $feat{$seq->id}{$primary_tag}{$locus_tag}=$feat;
    }
    
    # get meta data
    for my $feat(@feats){
      next if($feat->primary_tag()!~/source/i);
      $metaData{$seq->id}{source}=$feat;
      $metaData{$seq->id}{reference}=$reference;
    }
  }
  $in->close;
  
  return [\%feat,\%metaData];
}

sub readAnnotationDir{
  my($dir,$settings)=@_;
  my %fieldMap=annotationFieldMap();
  my %annotation;
  my $fileQueue=Thread::Queue->new(keys(%fieldMap));
  my @thr;
  local $$settings{numcpus}=1; # it's faster when it's not multithreaded for some reason. Disk I/O?
  $thr[$_]=threads->new(\&readSqlFile,$dir,$fileQueue,$settings) for(0..$$settings{numcpus}-1);
  $fileQueue->enqueue(undef) for(@thr);
  for(@thr){
    my $tmp=$_->join;
    while(my($locus_tag,$sqlAnnotations)=each(%$tmp)){
      while(my($sqlType,$annotationList)=each(%$sqlAnnotations)){
        $annotation{$locus_tag}{$sqlType}=$annotationList;
      }
    }
  }
  
  return \%annotation;
}


##################
## sub-subroutines
##################

sub readSqlFile{
  my($dir,$Q,$settings)=@_;
  my $tid="TID".threads->tid;
  my %fieldMap=annotationFieldMap();
  my %annotation;
  while(defined(my $mapKey=$Q->dequeue)){
    my @map=@{ $fieldMap{$mapKey} };
    next if(!grep/locus_tag/,@map); # don't process the sql if there is no locus_tag to match the annotation to
    # get the right filename for this sql file
    my $sqlfile="$dir/aa.fasta.$mapKey.sql";
    $sqlfile="$dir/$mapKey.sql" if(!-f $sqlfile);
    $sqlfile="$dir/aa.fasta.iprscan_out.xml.$mapKey.sql" if(!-f $sqlfile);
    logmsg "$tid $sqlfile";
    
    open(SQL,"sort $sqlfile | uniq |") or die "ERROR: could not open sql file $sqlfile: $!";
    while(<SQL>){
      next if(/^s*$/); # skip blank lines
      chomp;
      my %f;
      @f{@map}=readSqlLine($_);
      $f{evalue}=~s/,// if($f{evalue}); # get rid of those commas at the end (where do they come from??)
      push(@{ $annotation{$f{locus_tag}}{$mapKey} },\%f);
    }
    close SQL;
  }
  return \%annotation;
}

# TODO use uniprot and uniprot_evidence to add more useful information
sub interpretUniprot{
  my($geneFeat,$cdsFeat,$locusAnnotation,$settings)=@_;
  # take the best blast hit: highest score, lowest evalue
  $$locusAnnotation{blast}||=[];
  my @uniprotAnnotation=sort{
    return $$b{score}<=>$$a{score} if($$b{score}!=$$a{score});
    # not really sure why, but there is a problem with commas
    $$a{evalue}=~s/,//g;
    $$b{evalue}=~s/,//g;    
    return sprintf("%f",$$a{evalue})<=>sprintf("%f",$$b{evalue});
  } @{ $$locusAnnotation{blast} };
  my $uniprotAnnotation=$uniprotAnnotation[0] || {};
  my $uniprotDesc=$$uniprotAnnotation{description} || "";
  my $gene="";
  if($uniprotDesc=~/GN=([a-z]\S{2,7})/){ # gene name is 4+/-1 letters long. Probably never too much longer than that.
    $gene=$1;
  }
  my $proteinProduct="";
  if($uniprotDesc=~/^(^.+?(=|$))/){
    $proteinProduct=$1;
    $proteinProduct=~s/\S+=\S*$//;   # remove last word with equals sign
    $proteinProduct=~s/^\s+|\s+$//g; # trim whitespace
  }
  #my ($geneName,$proteinProduct,$uniprotDescription)=interpretUniprot($locusAnnotation,$settings);
  $geneFeat->add_tag_value('gene',$gene) if($gene);
  $cdsFeat->add_tag_value('product',$proteinProduct) if($proteinProduct);
  
  # add the uniprot evidence
  if($$uniprotAnnotation{hit_name}){
    my $note=($cdsFeat->get_tag_values('note'))[0];
    $cdsFeat->remove_tag('note');
    $note=~s/[\.;]?\s*$/. /;
    $note.="Product Predictor: Uniprot hit against $$uniprotAnnotation{hit_name}";
    $cdsFeat->add_tag_value('note',$note);
  
    # add the uniprot scores
    my $score="evalue: $$uniprotAnnotation{evalue}, bitscore: $$uniprotAnnotation{bits}";
    $cdsFeat->add_tag_value('score',$score);
  }
  
  return 1;
  
  # TODO think of something else if it hits against "putative", etc
  # Maybe go to the next hit, or use annotations from other tools.
  # TODO parse the description for more meaningful things
}

sub interpretIpr{
  my($cdsFeat,$iprAnnotation,$settings)=@_;
  my @newFeat;
  for my $an (@$iprAnnotation){
    # Figure out the correct start/stop
    my($ntCdsStart,$ntCdsEnd)=($$an{start}*3-3,$$an{end}*3-3); # aa to nt CDS coordinates, base 0
    my($start,$end)=($cdsFeat->start+$ntCdsStart, $cdsFeat->start+$ntCdsEnd); # CDS to genomic coordinates
    if($cdsFeat->strand<1){
      $end  =$cdsFeat->end-$ntCdsStart;
      $start=$cdsFeat->end-$ntCdsEnd;
    }
    
    my $newFeat=Bio::SeqFeature::Generic->new(-start=>$start,-end=>$end,
            -source_tag=>"InterPro",
            -primary=>"misc_feature",
            -strand=>$cdsFeat->strand,
            -tag=>{
              locus_tag=>($cdsFeat->get_tag_values('locus_tag'))[0],
            }
    );
    
    while(my ($key,$value)=each(%$an)){
      next if($key=~/^(start|end|locus_tag)$/); # exclude some tags from being shown like this
      $newFeat->add_tag_value($key,$value) if($value!~/^\s*$/); # who cares about blank values
    }
    
    $newFeat->start($start);
    $newFeat->end($end);
      
    push(@newFeat,$newFeat);
  }
  return @newFeat;
}

sub interpretSignalP{
  my($cdsFeat,$annotation,$settings)=@_;
  # currently: only paying attention to NNs and not HMMs
  my $signalp_nn=$$annotation{signalp_nn};
  my $is_signalpeptide=0;
  my ($start,$end,$score);
  for my $measure (@$signalp_nn){
    if($$measure{measure_type} eq 'D'){
      $is_signalpeptide=1 if($$measure{is_signal_peptide} eq 'YES');
      ($start,$end,$score)=($$measure{start},$$measure{end},$$measure{value});
    }
  }
  return 0 if(!$is_signalpeptide);
  
  my($ntCdsStart,$ntCdsEnd)=($start*3-3,$end*3-3); # 0-coordinates
  my($genomeStart,$genomeEnd)=($cdsFeat->start+$ntCdsStart,$cdsFeat->start+$ntCdsEnd);
  if($cdsFeat->strand<1){
    ($genomeEnd,$genomeStart)=($cdsFeat->end-$ntCdsStart,$cdsFeat->end-$ntCdsEnd);
  }
  my $signalpFeat=Bio::SeqFeature::Generic->new(-start=>$genomeStart,-end=>$genomeEnd,
            -source_tag=>"SignalP",
            -primary=>"sig_peptide",
            -strand=>$cdsFeat->strand,
            -tag=>{
              locus_tag=>($cdsFeat->get_tag_values('locus_tag'))[0],
              evidence=>"SignalP",
              score=>$score,
            }
  );
  return $signalpFeat;
}

sub interpretTmhmm{
  my ($cdsFeat,$annotation,$settings)=@_;
  my $tmhmm_location=$$annotation{tmhmm_location};
  my $tmhmm=$$annotation{tmhmm};
  $tmhmm=$$tmhmm[0]; # there's only one tmhmm result
  return 0 if(!$tmhmm);
  
  $tmhmm_location=[sort{$$a{start}<=>$$b{start}} @$tmhmm_location];
  my $splitLocation = new Bio::Location::Split();
  my (@start,@end);
  for(@$tmhmm_location){
    my $startNt=$$_{start}*3-3;
    my $endNt=$$_{end}*3-3;
    my $start=$cdsFeat->start+$startNt;
    my $end=$cdsFeat->start+$endNt;
    if($cdsFeat->strand<1){
      $start=$cdsFeat->end-$endNt;
      $end=$cdsFeat->end-$startNt;
    }
    $splitLocation->add_sub_Location(new Bio::Location::Simple(
        -start=>$start,
        -end=>$end,
        -strand=>$cdsFeat->strand,
    ));
  }
  
  my $tmFeat=Bio::SeqFeature::Generic->new(
    -primary=>"misc_structure",
    -source_tag=>"TMHMM",
    -strand=>$cdsFeat->strand,
    -location=>$splitLocation,
    -tag=>{
      locus_tag=>($cdsFeat->get_tag_values('locus_tag'))[0],
      evidence=>"TMHMM",
      score=>$$tmhmm{total_prob_n_in},
    }
  );
  
  $tmFeat->start(@start);
  $tmFeat->end(@end);
  return $tmFeat;
}

sub interpretIs{
  my($cdsFeat,$annotation,$settings)=@_;
  my $is=$$annotation{is_hits};
  return 0 if(!@$is);
  
  my $feat=blastSqlToFeat($cdsFeat,$is,"IS Finder database",{});
  return $feat;
}

sub interpretVf{
  my($cdsFeat,$annotation,$settings)=@_;
  my $vf=$$annotation{vfdb_hits};
  return 0 if(!@$vf);
  
  my $feat=blastSqlToFeat($cdsFeat,$vf,"VF database",{});
  return $feat;
}

sub interpretCogs{
  my($cdsFeat,$annotation,$extraAnnotationInfo,$settings)=@_;
  my $cogs=$$annotation{cogs_hits};
  return 0 if(!@$cogs);
  
  my $feat=blastSqlToFeat($cdsFeat,$cogs,"COGs database",{});
  my $cogsProt=($feat->get_tag_values('product'))[0];
  $feat->remove_tag('description');
  $feat->remove_tag('product');
  $feat->add_tag_value('description',$$extraAnnotationInfo{prot2cogs}{$cogsProt});
  return $feat;
}

##################
## tools
##################

# TODO do something with multiple results?
sub blastSqlToFeat{
  my($cdsFeat,$annotation,$evidence,$settings)=@_;
  my $a=$$annotation[0];
  my $feat=Bio::SeqFeature::Generic->new(
    -primary=>"misc_structure",
    -start=>$cdsFeat->start,
    -end=>$cdsFeat->end,
    -strand=>$cdsFeat->strand,
    -tag=>{
      evidence=>$evidence,
      score=>"evalue:$$a{evalue}, bitscore:$$a{bits}",
      locus_tag=>$$a{locus_tag},
      product=>$$a{hit_name},
      database=>$$a{db_name},
      description=>$$a{description},
    }
  );
}

# read a pipe-delimited sql line
sub readSqlLine{
  my($line,$settings)=@_;
  $line=~s/\\\|/::::/g; # protect all escaped pipes
  my @F=split(/\|/,$line);
  $_=~s/::::/|/g for(@F); # return all escaped pipes
  return @F;
}

sub uniprotSql{
  my($dir,$fieldMap,$settings)=@_;
  my @map=@{ $$fieldMap{uniprot} };
  my $uniprotFile="$dir/uniprot.sql";
  my %uniprot;
  open(UNIPROT,"sort $uniprotFile|uniq|") or die "Could not open uniprot evidence file at $uniprotFile: $!";
  while(<UNIPROT>){
    chomp;
    next if(/^s*$/); # skip blank lines
    my %f;
    @f{@map}=readSqlLine($_);
    $uniprot{$f{ac}}=\%f;
  }
  close UNIPROTEVIDENCE;
  return %uniprot;
}

sub uniprotEvidence{
  my($dir,$fieldMap,$settings)=@_;
  my $evFile="$dir/uniprot_evidence.sql";
  my @map=@{ $$fieldMap{uniprot_evidence} };
  my %evidence;
  open(UNIPROTEVIDENCE,"sort $evFile|uniq|") or die "Could not open uniprot evidence file at $evFile: $!";
  while(<UNIPROTEVIDENCE>){
    chomp;
    next if(/^s*$/); # skip blank lines
    my %f;
    @f{@map}=readSqlLine($_);
    $evidence{$f{ac}}=\%f;
  }
  close UNIPROTEVIDENCE;
  return %evidence;
}

# set up cogs mapping
sub prot2cogMapping{
  my($settings)=@_;
  my %prot2cogid;
  open(CFH,"<$$settings{prot2cogid}") or die "Could not open COGs mapping file $$settings{prot2cogid}: $!\nNeed prot2cogid in the conf file";
  while(<CFH>){
    chomp;
    my ($prot,$cogid)=split(/\s+/);
    $prot2cogid{$prot}=$cogid;
  }
  close CFH;
  return %prot2cogid;
}

# return a set of fields to understand annotation sql files
sub annotationFieldMap{
  #my $blastFields="locus_tag target_id evalue coverage db_name identity length description rank score bits percent_conserved hit_accession hit_name";
  my $blastFields="locus_tag target_id evalue coverage db_name identity length description rank score bits percent_conserved hit_accession hit_name";
  my @blastFields=split(/\s+/,$blastFields);
  my $blastFields2=$blastFields;
  $blastFields2=~s/target_id/uniprot_id/; #$blastFields2=~s/\bname\b/hit_name/;
  my @blastFields2=split(/\s+/,$blastFields2);
  my %map = (
    blast => [@blastFields2],
    vfdb_hits => [@blastFields],
    cogs_hits => [@blastFields],
    is_hits => [@blastFields],
    pdb_hits => [@blastFields],
    ipr_matches => [qw/locus_tag accession_num product database_name start end evalue status evidence/],
    signalp_hmm  => [ qw/locus_tag prediction signal_peptide_probability max_cleavage_site_probability start end/],
    signalp_nn => [ qw/locus_tag measure_type start end value cutoff is_signal_peptide/],
    tmhmm => [ qw/locus_tag length predicted_number expected_number_aa expected_number_aa_60 total_prob_n_in/],
    tmhmm_location => [qw/locus_tag location start end/],
    uniprot => [qw/ac ac2 source product pro_type gene_type gene_name gene_id/],
    uniprot_evidence => [qw/ac dbId db name/],
  );
  return %map;
}

sub usage{
  "Transforms a CG-Pipeline annotation directory into a standard genbank file
  Usage: ". basename($0) . " --prediction=prediction.gb --inputdir=annotation-sql-folder --gb=genbank-output-file --gff=gff-output-file [--organism=organism_id]
  "
}
