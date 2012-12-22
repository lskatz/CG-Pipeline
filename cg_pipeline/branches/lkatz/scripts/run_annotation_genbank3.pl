#!/usr/bin/env perl
use warnings;
use strict;
use FindBin;
use lib "$FindBin::RealBin/../lib";
use AKUtils qw/logmsg/;
use Bio::Seq;
use Bio::Seq::RichSeq;
use Bio::SeqIO;
use Bio::Tools::GFF;
use File::Basename;
use Data::Dumper;
use Getopt::Long;

exit(&main());

sub main(){
  my $settings={};
  die usage() if(@ARGV < 1);
  GetOptions($settings,qw(organism=s prediction=s inputdir=s gb=s gff=s)) or die;
  my $annotationDir=$$settings{inputdir};
  die "ERROR: cannot read annotation directory at $annotationDir" if(!-d $annotationDir || !-x $annotationDir);
  my $predictionGbk=$$settings{prediction};
  die "ERROR: cannot read prediction file at $predictionGbk" if(!-f $predictionGbk);
  my $outfile=$$settings{gb};
  die "ERROR: need genbank parameter".usage() if(!$outfile);
  
  logmsg "Transforming the prediction genbank into a set of features";
  my ($predfeatures,$metaData)=gbkToFeatures($predictionGbk,$settings);
  
  logmsg "Putting all annotations in $annotationDir into memory";
  my $annotation=readAnnotationDir($annotationDir,$settings);
  
  logmsg "Combining the features from the prediction file and the annotations";
  combinePredictionAndAnnotation($outfile,$predfeatures,$metaData,$annotation,$settings);
  
  return 0;
}

sub combinePredictionAndAnnotation{
  my ($outfile,$predfeatures,$metaData,$annotation,$settings)=@_;
  my %annotationDirMap=annotationFieldMap();
  my @annotationKeys=keys(%annotationDirMap);
  
  logmsg "One contig/chromosome per dot: ";
  my $out=Bio::SeqIO->new(-file=>">$outfile",-format=>"genbank");
  while(my ($seqid,$meta)=each(%$metaData)){
    $|++;print ".";$|--;
    my $sourceFeat=$$meta{source};
    my $seq=$sourceFeat->seq;
    $seq=Bio::Seq->new(-seq=>$seq->seq,-id=>$seq->id);
    $seq->add_SeqFeature($sourceFeat);
    my @sortedLocus=sort{
      my ($numA,$numB)=($a,$b);
      $numA=~s/\D+//g;
      $numB=~s/\D+//g;
      return $numA<=>$numB;
    } keys(%{$$predfeatures{$seqid}});
    
    # add annotations from the annotations directory
    for my $locus_tag(@sortedLocus){
      my $locusAnnotation=$$annotation{$locus_tag};
      if(!$locusAnnotation){
        $$locusAnnotation{$_}=[] for(@annotationKeys);
      }
      
      my $geneFeat=$$predfeatures{$seqid}{$locus_tag}->clone;
      my $cdsFeat =$$predfeatures{$seqid}{$locus_tag}->clone;
      $cdsFeat->primary_tag("CDS");
      my (
        @miscFeat,  # for signalP, tmhmm, and other "named" tags in the gb file
        @otherFeat, # for protein domains
      );
      
      # modify the gene and cds features by reference
      interpretUniprot($geneFeat,$cdsFeat,$locusAnnotation,$settings);
      
      # signal peptide
      my $signalpFeat=interpretSignalP($cdsFeat,$locusAnnotation,$settings);
      push(@miscFeat,$signalpFeat) if($signalpFeat); # $signalpFeat is 0 if not present
      
      # transmembrane helix
      my $tmFeat=interpretTmhmm($cdsFeat,$locusAnnotation,$settings);
      push(@otherFeat,$tmFeat) if($tmFeat);
      
      # COGs
      
      # others: IS elements, VFs, ...
      
      # domain hits
      my @iprFeat=interpretIpr($cdsFeat,$$locusAnnotation{ipr_matches},$settings);
      push(@otherFeat,@iprFeat);
      
      @otherFeat=sort {$a->start<=>$b->start} @otherFeat;
      $seq->add_SeqFeature($geneFeat,$cdsFeat,@miscFeat,@otherFeat);
    }
    
    $out->write_seq($seq);
  }
  $out->close;
  print "\n"; # newline after all the "update dots"
  
  return 1;
}

sub gbkToFeatures{
  my($gbk,$settings)=@_;
  my %feat;
  my %metaData;
  my $in=Bio::SeqIO->new(-file=>$gbk);
  while(my $seq=$in->next_seq){
    my @feats = $seq->all_SeqFeatures();
    
    # looking at gene features only
    for my $feat(@feats){
      #delete($$feat{'_gsf_seq'}); # save memory and simplify things by removing the embedded seq object (hopefully not buggy)
      next if($feat->primary_tag()!~/gene/i); # only look at genes
      next if(!$feat->has_tag('locus_tag'));
      my @locus_tag=$feat->get_tag_values('locus_tag');
      die "ERROR: a locus in $gbk has more than one locus_tag" if(@locus_tag>1);
      my $locus_tag=$locus_tag[0];
      die "ERROR: locus tag $locus_tag exists twice." if($feat{$seq->id}{$locus_tag});
      $feat{$seq->id}{$locus_tag}=$feat;
    }
    # get meta data
    # TODO add species info; author info; etc
    for my $feat(@feats){
      next if($feat->primary_tag()!~/source/i);
      $metaData{$seq->id}{source}=$feat;
    }
  }
  $in->close;
  return (\%feat,\%metaData);
}

# TODO multithread this
sub readAnnotationDir{
  my($dir,$settings)=@_;
  my %fieldMap=annotationFieldMap();
  my %uniprotEvidence=uniprotEvidence($dir,\%fieldMap,$settings);
  my %uniprot=uniprotSql($dir,\%fieldMap,$settings);
  my @sqlfile=map("aa.fasta.$_.sql",keys(%fieldMap));
  my %annotation;
  for my $mapKey(keys(%fieldMap)){
    my @map=@{ $fieldMap{$mapKey} };
    next if(!grep/locus_tag/,@map); # don't process the sql if there is no locus_tag to match the annotation to
    # get the right filename for this sql file
    my $sqlfile="$dir/aa.fasta.$mapKey.sql";
    $sqlfile="$dir/$mapKey.sql" if(!-f $sqlfile);
    $sqlfile="$dir/aa.fasta.iprscan_out.xml.$mapKey.sql" if(!-f $sqlfile);
    logmsg $sqlfile;
    
    open(SQL,"sort $sqlfile | uniq |") or die "ERROR: could not open sql file $sqlfile: $!";
    while(<SQL>){
      next if(/^s*$/); # skip blank lines
      chomp;
      my %f;
      @f{@map}=readSqlLine($_);
      push(@{ $annotation{$f{locus_tag}}{$mapKey} },\%f);
    }
    close SQL;
  }
  return \%annotation;
}


##################
## sub-subroutines
##################

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
  return 1;
  
  # TODO think of something else if it hits against "putative", etc
  # Maybe go to the next hit, or use annotations from other tools.
  # TODO parse the description for more meaningful things
}

sub interpretIpr{
  my($cdsFeat,$iprAnnotation,$settings)=@_;
  my @newFeat;
  for my $an (@$iprAnnotation){
    my $newFeat=$cdsFeat->clone;
    $newFeat->primary_tag("misc_feature");
    while(my ($key,$value)=each(%$an)){
      next if($key=~/^(start|end|locus_tag)$/); # exclude some tags from being shown like this
      $newFeat->add_tag_value($key,$value) if($value!~/^\s*$/); # who cares about blank values
    }
    
    # Figure out the correct start/stop
    my($ntCdsStart,$ntCdsEnd)=($$an{start}*3-3,$$an{end}*3-3); # aa to nt CDS coordinates, base 0
    my($start,$end)=($cdsFeat->start+$ntCdsStart, $cdsFeat->start+$ntCdsEnd); # CDS to genomic coordinates
    if($newFeat->strand<1){
      $end  =$cdsFeat->end-$ntCdsStart;
      $start=$cdsFeat->end-$ntCdsEnd;
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

##################
## tools
##################
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
