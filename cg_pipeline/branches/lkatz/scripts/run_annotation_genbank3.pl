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
  
  my $out=Bio::SeqIO->new(-file=>">$outfile",-format=>"genbank");
  while(my ($seqid,$meta)=each(%$metaData)){
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
      my @otherFeat;
      
      # modify the gene and cds features by reference
      interpretUniprot($geneFeat,$cdsFeat,$locusAnnotation,$settings);
      
      # domain hits
      my @iprFeat=interpretIpr($cdsFeat,$$locusAnnotation{ipr_matches},$settings);
      push(@otherFeat,@iprFeat);
      
      # signal peptide
      my ($is_signalPeptide,$spStart,$spEnd)=interpretSignalP($locusAnnotation,$settings);
      $spStart+=$geneFeat->start-1;
      $spEnd  +=$geneFeat->start-1;
      # TODO what to do for revcomp sequences???
      if($is_signalPeptide){
        my $signalpFeat=Bio::SeqFeature::Generic->new(-start=>$spStart,-end=>$spEnd,
            # TODO: -score, 
            -source_tag=>"SignalP",
            -primary=>"sig_peptide",
            -strand=>$geneFeat->strand,
            -tag=>{
              locus_tag=>$locus_tag,
            }
        ); 
        push(@otherFeat,$signalpFeat);
      }
      
      # transmembrane helix
      
      # COGs
      
      # others: IS elements, VFs, ...
      
      $seq->add_SeqFeature($geneFeat,$cdsFeat,@otherFeat);
    }
    
    $out->write_seq($seq);
  }
  
  $out->close;
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
  if($uniprotDesc=~/GN=([a-z]\S{2,9})/){ # gene name is 3+/-1 letters long. Probably never too much longer than that.
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
    $newFeat->primary_tag("misc_structure");
    while(my ($key,$value)=each(%$an)){
      next if($key=~/^(start|end|locus_tag)$/); # exclude some tags from being shown like this
      $newFeat->add_tag_value($key,$value) if($value!~/^\s*$/); # who cares about blank values
    }
    
    # TODO I don't think that the start/stop sites are exactly on the domain.
    # I think that they are start/stop for the gene which is wrong
    my($start,$end)=($cdsFeat->start+$$an{start}-1, $cdsFeat->end+$$an{start}-1);
    if($newFeat->strand<1){
      $start=$cdsFeat->end-($end-$start);
      $end=$cdsFeat->end;
    }
    $newFeat->start($start);
    $newFeat->end($end);
      
    push(@newFeat,$newFeat);
  }
  return @newFeat;
}

sub interpretSignalP{
  my($annotation,$settings)=@_;
  # currently: only paying attention to NNs and not HMMs
  my $signalp_nn=$$annotation{signalp_nn};
  my $is_signalpeptide=0;
  my ($start,$end);
  for my $measure (@$signalp_nn){
    if($$measure{measure_type} eq 'D'){
      $is_signalpeptide=1 if($$measure{is_signal_peptide} eq 'YES');
      ($start,$end)=($$measure{start},$$measure{end});
    }
  }
  return ($is_signalpeptide,$start,$end);
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
