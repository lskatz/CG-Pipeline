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
  if(@ARGV < 1){
    print STDERR "Usage: ". basename($0) . " --prediction=prediction.gb --inputdir=annotation-sql-folder --gb=genbank-output-file --gff=gff-output-file [--organism=organism_id] \n";
    return 0;
  }
  my %data;
  my @args;
  my $newftr;
  my $tag;
  my $settings = {
    appname => 'cgpipeline',
  };
  # blastFields is for any output from run_annotation_blast.pl
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
  );

  $settings = AKUtils::loadConfig($settings);
  GetOptions($settings,('organism=s','prediction=s','inputdir=s','gb=s','gff=s')) or die;

  # COGs mapping
  my %prot2cogid=prot2cogMapping($settings);
  my $organism=(defined($$settings{'organism'}))?$$settings{'organism'}:"organism";
  
  # load up the correct SQL files from the annotation directory
  my @sqlfiles;
  my $atndir;
  opendir($atndir,$$settings{'inputdir'}) or die "unable to open directory $$settings{'inputdir'}:$!\n";
  my @atnfiles=readdir($atndir);
  foreach my $type ((qw/blast ipr_matches signalp_hmm signalp_nn tmhmm.sql tmhmm_location vfdb_hits cogs_hits/)){ # I removed uniprot so that it could be specially parsed
    my @files=grep(/$type/,@atnfiles);
    if(@files){push(@sqlfiles,$$settings{'inputdir'} . "/" . $files[0]);}
  }
  
  # Set up I/O
  my $gbin = Bio::SeqIO->new(-file=>$$settings{'prediction'},-format=>'genbank');
  my $gbout = Bio::SeqIO->new(-file=>">".$$settings{'gb'},-format=>'genbank');
  my $fasta;
  if(defined($$settings{'fasta'})){
    $fasta = Bio::SeqIO->new(-file=>">".$$settings{'fasta'},-format=>'');
  }
  my $gff;
  if(defined($$settings{'gff'})){
    $gff = Bio::Tools::GFF->new(-file=>">".$$settings{'gff'},-gff_version => 3);
  }

  STAGE_1:
  #Stage 1: Make a hash of loci to SeqFeature objects from each data line
  my %uniprot=uniprotInfo("$$settings{inputdir}/uniprot.sql",$map{uniprot},$settings);
  my %spfeats=annotationsHash(\@sqlfiles,\%map,$settings);

  STAGE_2:
  logmsg "Stage 2: Select a name for each locus from available SeqFeature objects and integrate the other subfeatures";
  while( my $seq=$gbin->next_seq()){ #each contig
    logmsg "Working on ".$seq->id;
    if($organism eq 'organism'){#not set
      $organism=$seq->primary_seq->desc;
      $organism =~ s/^([^,_]+).*$/$1/;
    }
    my $contig=$seq->primary_seq->display_id;
    $contig =~ s/\D*[a-zA-Z]*([0-9]+).*$/$1/;
    $seq->display_id(sprintf("%s_%04d",$organism,$contig));
    if(defined($$settings{'pipeline_version'})){#write pipeline version into the COMMENT field
      my $pipeline_version="CG-Pipeline version " . $$settings{'pipeline_version'};
      $seq->annotation->add_Annotation('comment',Bio::Annotation::Comment->new(-text=>$pipeline_version));
    }
    if(defined($$settings{'pipeline_reference'})){
      my ($gb_authors,$gb_title,$gb_journal,$gb_pubmed)=split('\|',$$settings{'pipeline_reference'});
      $seq->annotation->add_Annotation('reference',Bio::Annotation::Reference->new(-title=>$gb_title,-authors=>$gb_authors,-location=>$gb_journal,-pubmed=>$gb_pubmed));
    }
    my @feats = $seq->all_SeqFeatures(); # features here are loci in the genome
    #$seq->flush_SeqFeatures(); # try removing this.....
    # default source tag if it's not in the prediction file
    my $sourceTag=Bio::SeqFeature::Generic->new(
      -start=>1,
      -end=>$seq->length,
      -primary=>"source",
      -tag=>{
        organism=>"organismXYZ",
      },
    );

    # Find all annotations and add them to this feature.
    # Add the new feature to the output sequence
    for my $predFeat(@feats){
      my $locus_tag = ($predFeat->get_tag_values('locus_tag'))[0];
      
      my @subfeats; # all features for this gene
    }
    
    $gbout->write_seq($seq);
  }
  
  logmsg "Done!";
  return 0;
}#end main

sub uniprotInfo{
  my($uniprot,$uniprotMap,$settings)=@_;
  my %uniprot;
  open (FH, "<$uniprot") or die "Could not open $uniprot: $!";
  while(<FH>){
    chomp;
    # skip a row of empties
    if( /^\s*\|/ ){ print STDERR "Skipping blank line:$_\n"; next;}
    s/\\\|/::/g; # escaped pipes are replaced
    my @values = split(/\|/);
    @{$uniprot{$values[0]}}{@$uniprotMap}=@values;
  }
  close FH;
  return %uniprot;
}

sub annotationsHash{
  my($sqlfiles,$map,$settings)=@_;
  my %map=%$map;
  my $blastcount=0;
  my %spfeats;

  # get annotations from all annotation sources
  my %cgptype;
  foreach my $sqlfile(@$sqlfiles){
    my($sqlname,$sqlpath,$sqlsuffix)=fileparse($sqlfile,qw(.sql));
    my $type=$sqlname; $type=~s/^.+\.//;
    if(!defined($type) || !defined($map{$type})){ print STDERR "Input filename rejected: $type\n";next;}
    $cgptype{basename($seqlfile)}=1;
    my @params = @{$map{$type}};
    open (FH, "<$sqlfile") or die $!; #one prediction file containing locus tags
    logmsg "Parsing $sqlfile";
    while (<FH>){
      chomp;
      # skip a row of empties
      if( /^\s*\|/ ){ print STDERR "Skipping blank line:$_\n"; next;}
        
      my $newftr = Bio::SeqFeature::Generic->new(-primary=>$type,-start=>'0',-end=>'0');
      $newftr->add_tag_value("cgptype",$type);
      s/\\\|/::/g;
      my @values = split(/\|/);
      if(scalar @values != scalar @params){print STDERR "Wrong number of fields in $sqlfile:" . scalar @values . "\n"."  ".join("____",@values)."\n";next;}

      
      for(my $i=0;$i<scalar @params;$i++){
        if($params[$i] eq "start"){
          $newftr->start($values[$i]);
        }
        elsif($params[$i] eq "end"){
          $newftr->end($values[$i]);
        }
        else{
          $newftr->add_tag_value($params[$i],$values[$i]);
        }
      }
      if(!defined($spfeats{$values[0]})){$spfeats{$values[0]} = [];}
      push (@{$spfeats{$values[0]}}, $newftr );
    }
    close(FH);
  }

  # now add in the information given from prediction.gb
  my $in=Bio::SeqIO->new(-file=>$$settings{'prediction'},-format=>"genbank");
  while(my $seq=$in->next_seq){
    my @feats = $seq->all_SeqFeatures();
    for my $ftr(@feats){
      next if(!$ftr->has_tag('locus_tag')); # source tag, etc
      my $locus_tag = ($ftr->get_tag_values('locus_tag'))[0];
      push(@{$spfeats{$locus_tag}},$ftr);
    }
  }

  # TODO load in cgptype?
  %spfeats=combineSpFeats(\%spfeats,$settings);

  return %spfeats;
}

sub combineSpFeats{
  my($spfeats,$settings)=@_;
  my %newFeats;
  while(my($featid,$featArr)=each(%$spfeats)){
    my $mergedFeatArr=combineSpFeat($featArr,$settings);
    $newFeats{$featid}=$mergedFeatArr;
  }
  return %newFeats;
}
sub combineSpFeat{
  my($featArr,$settings)=@_;
  print Dumper $featArr;exit;
  my %f; # feature hash with keys=>gene, CDS, etc
  for my $f(@$featArr){
    my $cgptype=($f->get_tag_values('cgptype'))[0];
    die "TODO reconstruct this feature if it doesn't already exist in %f";
  }
}

# set up cogs mapping
sub prot2cogMapping{
  my($settings)=@_;
  my %prot2cogid;
  open(CFH,"<$$settings{prot2cogid}") or die;
  my @prot2cogid_map=<CFH>;
  close (CFH);
  foreach (@prot2cogid_map){
    my ($prot,$cogid)=split(/\s+/);
    $prot2cogid{$prot}=$cogid;
  }
  undef @prot2cogid_map;
  return %prot2cogid;
}

sub gene_factory($$$$){
  my ($start,$end,$strand,$tags) = @_;
  return Bio::SeqFeature::Generic->new(-primary=>'gene',-source_tag=>'UniProt',-start=>$start,-end=>$end,-strand=>$strand,-tag=>$tags);
}
sub remove_tags($){
  my $ftr = shift;
  my @note=();
  my @tags=$ftr->get_all_tags();
  my $value;
  foreach my $tag (@tags){
    $value=($ftr->get_tag_values($tag))[0];
    $value=~s/=/:/g;
    @note = (@note,
        ("$tag:" . $value)
      );
  }
  foreach my $tag (@tags){
    $ftr->remove_tag($tag);
  }
  return join("; ",@note);
}
sub replace_tags($){
  my $ftr = shift;
  if(!$ftr->has_tag('note')){return;}
  my @note = split(';',join('; ',$ftr->get_tag_values('note')));
  foreach (@note){
    my ($tag,$value)=split(':');
    if($ftr->has_tag($tag)){next;}#things like 'product'
    else{$ftr->add_tag_value($tag,$value);}
  }
  $ftr->remove_tag('note');
  my @tags = $ftr->all_tags();
  #fix bad tag names
  foreach my $tag (@tags){
    my @values = $ftr->get_tag_values($tag);
    $ftr->remove_tag($tag);
    if(0<scalar @values){next;}
    $tag =~ s/^\s+//;
    #remove duplicate values
    my @uniquevalues=("uninitialized");
    foreach my $value(sort @values){
      if($value ne $uniquevalues[$#uniquevalues]){push(@uniquevalues,$value);}
      else{push(@uniquevalues,$value);}
    }  
    if(0<scalar @uniquevalues){
      $ftr->add_tag_value($tag,@uniquevalues);
    }
  }
}
sub start2nuc($$){
  my ($start,$coord) = @_;
  return $start - 3 + 3*$coord;
}
sub end2nuc($$){
  my ($start,$coord) = @_;
  return $start - 1 + 3*$coord;
}
sub tags2hash($){
  my $ftr = shift;
  my %data = ();
  my @tags=$ftr->get_all_tags();
  foreach my $tag (@tags){$data{$tag}=($ftr->get_tag_values($tag))[0];}
  return \%data;
}
sub addtranslation($$){
  my $ftr = shift;
  my $seq = shift;
  if($ftr->strand < 0){
    $ftr->add_tag_value('translation',Bio::Seq->new(-seq=>$seq->subseq($ftr->start,$ftr->end))->revcom->translate->seq);
  }
  else{
    $ftr->add_tag_value('translation',Bio::Seq->new(-seq=>$seq->subseq($ftr->start,$ftr->end))->translate->seq);
  }
}
