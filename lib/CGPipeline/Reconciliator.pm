#!/usr/bin/env perl

# CGPipeline::Reconciliator: Run Zimin et al's Reconciliator, create proper input files, handle output files
# Author: Lee Katz <lkatz@cdc.gov>

package CGPipeline::Reconciliator;
require 5.005;
my ($VERSION) = ('$Id$' =~ /,v\s+(\d+\S+)/o);

use strict;
use Exporter;
use List::Util qw(min max sum reduce shuffle);
use Data::Dumper;
use Bio::Perl;
use Bio::Assembly::IO;
use File::Path;
use File::Basename;
use File::Temp /qw tempfile tempdir /;
use File::Spec;
use threads;
use Thread::Queue;

local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };

our @ISA = "Exporter";
our @methods = qw();
our %EXPORT_TAGS = (all => [@methods]);
Exporter::export_ok_tags('all');

sub logmsg {my $FH = $AKUtils::LOG || *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

sub new {
	my ($class,$args) = @_;
	$class = ref($class) || $class;
  my $self=$args;
  $$self{inputDir}=[] if(!$$self{inputDir});
  $$self{num_bases_per_alignment_batch}||=50000000;
  $$self{cpus}||=1;
  $$self{default_good_qual_value}||=70;
  bless $self, $class;
  return $self;
}

sub reconcile{
  my($self,$args)=@_;
  
  my $executable=AKUtils::fullPathToExec("reconcile2.sh");
  my $executableDir=dirname($executable);
  my @rDir=@{ $self->{inputDir} };
  my $numDirs=@rDir;
  for(my $i=0;$i<$numDirs;$i++){$rDir[$i]=File::Spec->rel2abs($rDir[$i]);}
  $numDirs=2; logmsg "Number of assemblies has been artificially set to 2";
  # TODO sort by fewest contigs to get the best reference assembly and to put the
  # supplementary assemblies last
  for(my $i=0;$i<$numDirs;$i++){
    my $referenceAssembly=$rDir[$i];
    my $mateinfoFile="$referenceAssembly/mateinfo.txt";
    for(my $j=$i+1;$j<$numDirs;$j++){
      my $supplementaryAssembly=$rDir[$j];
      my $command="cd $executableDir; ";
      $command.="$executable $referenceAssembly $supplementaryAssembly $mateinfoFile ".$self->{num_bases_per_alignment_batch}." ".$self->{cpus};
      $command.=" 2>&1";
      logmsg "COMMAND $command";
      system($command);
    }
  }
  die;
}

# in order for the Reconciliator to work it must have the following files
# contigs.bases - contigs fasta file
# contigs.quals - contigs quality file (fasta format with phred scores)
# supercontigs
#   supercontig s1
#   contig c1
#   gap 200 * * 2
#   contig c7
#   gap 2235 * * 5
#   contig c22
# reads.placed - has the following fields
#   NCBI TI number (just put * for unknown)
#   read name
#   start of trimmed read (1-based coordinates in original-read-space)
#   number of bases in trimmed read
#   orientation on contig (0=forward, 1=reverse)
#   contig name
#   supercontig name
#   approximate start of trimmed read on contig (1-based coordinates)
#   approximate start of trimmed read on supercontig (1-based coordinates)
# reads.unplaced - has two parts of the same file
#   first part: short explanation mapped to a description of the explanation
#     chimera: "suspected of being chimeric"
#     repeat: "Probably part of a repeat element"
#   second part: TI number to short reason mapping
#     TI343343 Repeat
#     TI489568 Chimera

############################
###### raw Reconciliator data
############################
sub addReconciliatorDirectory{
  my($self,$args)=@_;
  my $rDir=$$args{rDir} or die "Error: 'rDir' parameter is needed.";
  
  # check required files
  for my $file(qw(reads.placed reads.unplaced contigs.quals contigs.bases)){
    if(!-f "$rDir/$file"){
      warn "Could not add directory $rDir because it does not exist!\n";
      return 0;
    }
  }

  # check optional files
  for my $file (qw(supercontigs mateinfo.txt)){
    if(!-f "$rDir/$file"){
      system("touch $rDir/$file");
    }
  }

  push(@ {$self->{inputDir} },$rDir);
  return 1;
}

############################
###### Velvet
############################
sub addVelvetAssembly($$){
  my($self,$args)=@_;
  my $velvetBaseDir=$$args{velvetDir} || die "Need 'velvetDir' argument.";
  my $baseDir="$$args{velvetDir}/reconciliator";
  $$args{baseDir}=$baseDir;
  $$args{readsfasta}="$$args{velvetDir}/Sequences";
  system("mkdir $baseDir") if(!-e $baseDir);
  return $self->_convertVelvetToSangerWashu($args);
}
sub _convertVelvetToSangerWashu($$){
  my($self,$args)=@_;

  my $contigsFile="$$args{velvetBaseDir}/contigs.fa";
  my $qualFile=""; # TODO create a fake quality file "$velvetBaseDir/contigs.qual";
  my $supercontigsFile=""; # TODO if exists
  $$args{acePath}=$self->_afgToAce("$$args{velvetDir}/velvet_asm.afg");
  $self->_convertAce($args);
  die "Need to make libraries.txt and read.pairs if possible: http://www.genome.umd.edu/reconciliation_instructions.htm. Probably make the libraries.txt file from the velvet afg parser in velvet/contrib folder. Read.pairs can probably be made by something in the afg/ace file";
}

sub _afgToAce{
  my($self,$afg)=@_;
  my ($filename,$dir)=fileparse($afg);
  my $ace=$afg; $ace=~s/\.afg$/.ace/;
  if(-e $ace){
    logmsg "$ace already exists. Skipping its creation. `rm $ace` to force its creation next time.";
  } else {
    logmsg "Converting $afg to $ace";
    system("amos2ace $afg");
  }
  return $ace;
}

##################################
#### Generic
##################################

# see if a contig is probably a velvet one (i.e. no meaningful qual scores)
sub _isVelvet{
  my($self,$contig)=@_;
  my $is_velvet=0;
  my @velvetQualScores=qw(17 19 23 36);
  #if the quality scores only have 17, 19, 23, and 36, then it is velvet
  my $qualityScoreStr=$contig->get_consensus_quality->seq;
  my @quality=split(/\s+/,$qualityScoreStr);
  @quality=splice(@quality,0,1000); # sample the first thousand
  my %seen=();
  for my $qualValue(@quality){
    $seen{$qualValue}=1;
  }
  for my $v (@velvetQualScores){
    if(!$seen{$v}){ last;}
    delete($seen{$v});
  }
  # it's a velvet assembly if the key qual values were removed and nothing remains
  if(scalar(keys(%seen))<1){
    $is_velvet=1;
  }
  return $is_velvet;
}

sub _removeGapsFromContig{
  my($self,$contig,$args)=@_;
  my $gap='-';
  my $Seq=$contig->get_consensus_sequence;
  my $Qual=$contig->get_consensus_quality;

  my $newSeq=Bio::Seq->new(-id=>$Seq->id,-desc=>$Seq->desc, -seq=>$Seq->seq);
  my $newQual=Bio::Seq::PrimaryQual->new(-id=>$Qual->id,-desc=>$Qual->desc, -qual=>$Qual->seq);

  # change each qual value to $$args{default_good_qual_value} if velvet assembly
  if($$args{is_velvet}){
    my $qualStr=$newQual->seq;
    my $newQualValue=$self->{default_good_qual_value};
    $qualStr=~s/\d+/$newQualValue/g;
    $newQual->seq($qualStr);
  }
  return ($newSeq,$newQual);
}

# requires latest bioperl which is why the lib includes the latest Bio::Assembly (2010-10-24)
# return 0 or the directory where the conversion is stored
sub _convertAce($$){
  my($self,$args)=@_;
  logmsg "Converting $$args{acePath} into Reconciliator format in $$args{baseDir}";
  my $readsPlacedFilename="$$args{baseDir}/reads.placed";
  my $readsUnplacedFilename="$$args{baseDir}/reads.unplaced";
  my $contigsBases=Bio::SeqIO->new(-format=>"fasta",-file=>">$$args{baseDir}/contigs.bases");
  my $contigsQuals=Bio::SeqIO->new(-format=>"qual",-file=>">$$args{baseDir}/contigs.quals");

  # create a mapping hash from velvet read names to original read names
  # TODO if not velvet, don't do this mapping
  # TODO detect if it is velvet, if the first contig only has quality scores of "17 19 23 36", or put it in settings somehow
  #my $readMap={};
  #my $reverseReadMap={};
  #my $command="grep '>' $$args{velvetDir}/Sequences|cut -f1,2"; # returns the read original read name and the read name of each sequence
  #my $out=`$command`;
  #my @line=split(/\n/,$out);
  #for(my $i=0;$i<@line;$i++){
  #  if($line[$i]=~/>?(\S+).+?(\d+)$/){
  #    my($origName,$newName)=($1,$2);
  #    $$readMap{$newName}=$origName;
  #    $$reverseReadMap{$origName}=$newName;
  #  }
  #}
  #print Dumper $readMap;exit;

  # make a hash of all read names
  my $readName={};
  logmsg "Reading reads file $$args{readsfasta}";
  my $seqio=Bio::SeqIO->new(-file=>$$args{readsfasta},-format=>"fasta");
  while(my $seq=$seqio->next_seq){
    $$readName{$seq->id}=1;
  }
  my $ace=$$args{acePath};
  my $aceio=Bio::Assembly::IO->new(-file=>$ace,-format=>"ace");
  logmsg "Making reads placed and reads unplaced files";
  open(READSPLACED,">$readsPlacedFilename") or die("Could not open $readsPlacedFilename for writing");
  open(READSUNPLACED,">$readsUnplacedFilename") or die("Could not open $readsUnplacedFilename for writing");
  #TODO multithread based on each contig
  while(my $contig=$aceio->next_contig){
    logmsg "Converting contig ". $contig->id;
    
    $$args{is_velvet}=$self->_isVelvet($contig) if(!defined($$args{is_velvet}));

    my ($contigSeq,$contigQual)=$self->_removeGapsFromContig($contig,$args);
    $contigsBases->write_seq($contig->get_consensus_sequence);
    $contigsQuals->write_seq($contig->get_consensus_quality);
    warn "Having problems here in making a contigs file without gaps";
    if(ref($contig) eq 'Bio::Assembly::Contig'){
      my @read=$contig->get_seq_ids; 
      for my $r (@read){
        $r=$contig->get_seq_by_name($r);
        my $readId=$r->id;
        delete($$readName{$readId});
        my $numBases=abs($r->end-$r->start)+1;
        my $supercontigName='*';
        my $strand=($r->strand>0)?1:0;
        $strand='*' if($r->strand==0);
        my $contigStart=min($r->start,$r->end);
        my $contigStop=max($r->start,$r->end);
        my $readStart=1; #this should be where the read was trimmed but I cannot figure this out yet
        my $readStop=$r->length;
        print READSPLACED join(" ",$readId,$readId,$readStart,$readStop,$numBases,$strand,$contig->id,$supercontigName,$contigStart,$contigStop);
        print READSPLACED "\n";
      }
    }
    elsif(ref($contig) eq 'Bio::Assembly::Singlet'){
      my @read=$contig->get_seq_ids;
      my $r=$contig->get_seq_by_name($read[0]);
      my $readId=$r->id;
      delete($$readName{$readId});
      print READSUNPLACED join(" ",$readId,$readId,"singlet")."\n";
    }
    else{
      logmsg "Warning: Encountered new contig type, ".ref($contig).". Unsure of what to do. Skipping.";
    }
  }

  # get remaining reads and put them into unplaced.reads
  # TODO get the real reason why it wasn't incorporated, via the ace file
  foreach my $readId (keys %$readName){
    print READSUNPLACED join(" ",$readId,$readId,'singlet')."\n";
  }
  close READSPLACED;
  close READSUNPLACED;
  $self->addReconciliatorDirectory({rDir=>$$args{baseDir}});
  return $$args{baseDir};
}

############################
#### generic ace file
############################
sub addAce($$){
  my($self,$args)=@_;
  $$args{acePath}=$$args{ace} || die "Need 'ace' argument.";
  my $acedir=dirname($$args{acePath});
  $$args{baseDir}=File::Temp::tempdir("$$args{tempdir}/Reconciliator_generic_XXXXX",CLEANUP=>0, DIR=>$acedir);
  
  if(! -e "$$args{baseDir}/reads.placed" || $$args{force}){
    logmsg "Converting $$args{acePath} to Reconciliator input";
    system("mkdir $$args{baseDir}") if(!-e $$args{baseDir});
    $self->_convertAce($args);
  } else {
    logmsg "The Reconciliator input files already exist--skipping the conversion. To force conversion, use --force";
  }
  return $$args{baseDir};
}

############################
#### newbler
############################
sub addNewblerAce($$){
  my($self,$args)=@_;
  my $newblerAce=$$args{newblerAce} || die "Need 'newblerAce' argument.";
  my $acedir=dirname($newblerAce);
  $$args{baseDir}=File::Temp::tempdir("$$args{tempdir}/Reconciliator_newbler_XXXXX",CLEANUP=>0, DIR=>$acedir);
  $$args{ace454Path}=$newblerAce;
  $$args{acePath}="$$args{baseDir}/Contigs.ace";
  logmsg "Converting $$args{ace454Path} to a standard ace file.";
  $self->_newblerAceToAce($args);

  if(! -e "$$args{baseDir}/reads.placed" || $$args{force}){
    logmsg "Converting $$args{acePath} to Reconciliator input";
    system("mkdir $$args{baseDir}") if(!-e $$args{baseDir});
    $self->_convertAce($args);
  }
  else{
    logmsg "The Reconciliator input files already exist--skipping the conversion. To force conversion, use --force";
  }
  return $$args{baseDir};
}
sub _newblerAceToAce($args){
  my($self,$args)=@_;
  my $ace454=Bio::Assembly::IO->new(-file=>$$args{ace454Path},-format=>"ace",-variant=>'454');
  my $ace=Bio::Assembly::IO->new(-file=>">$$args{acePath}",-format=>"ace"); #output ace
  my $numContigs=`grep -c ^CO $$args{ace454Path}`+0;
  logmsg "Converting $$args{ace454Path} (454-ace) to $$args{acePath} (ace). $numContigs contigs total.";
  $|++; # get the progress bar in real-time by adding one to caching
  while(my $contig=$ace454->next_contig){
    print "."; # just print out a dot per contig as a progress meter
    $ace->write_contig($contig);
  }
  print "\n"; # newline after all those dots
  $|--; # put caching back where it was
  return $$args{acePath};
}

sub _createReadsUnplacedFile($$){
  my($self,$args)=@_;
  my($readHash);
  my $outputFilename="$$args{baseDir}/reads.unplaced";
  open(OUT,">",$outputFilename) or die("Could not open $outputFilename for writing");
  #   first part: short explanation mapped to a description of the explanation
  #     chimera: "suspected of being chimeric"
  #     repeat: "Probably part of a repeat element"
  my %explanationMapping=(
    Assembled=>"The read is fully incorporated into the assembly.",
    PartiallyAssembled=>"Only part of the read was included in the assembly; the rest was deemed to have divered sufficiently to not be included.",
    Singleton=>"The read did not overlap with any other reads in the input.",
    Repeat=>"The read was identified by the assembler as likely coming from a repeat region, and so was excluded from the final contigs.",
    Outlier=>"The read was identifi ed by the GS De Novo Assembler as problematic, and was excluded from the fi nal contigs (one explanation of these outliers are chimeric sequences, but sequences may be identifi ed as outliers simply as an assembler artifact).",
    TooShort=>"the trimmed read was too short to be used in the computation (i.e., shorter than 50 bases, unless 454 Paired End Reads are included in the dataset, in which case, shorter than 15 bases).",
  );
  while( my($k,$v)=each(%explanationMapping)){
    print OUT "$k:\"$v\"\n";
  }

  #   second part: TI number to short reason mapping
  #     TI343343 Repeat
  #     TI489568 Chimera
  my $read454Reads=$self->_read454ReadStatusFile("$$args{newblerDir}/454ReadStatus.txt");
  while( my($k,$v)=each(%$read454Reads)){
    next if($$v{status}=~/Assembled|Full|Partial/i);

    print OUT join("\t",($$v{TI},$$v{status}))."\n";
  }
  close OUT;
  return $outputFilename;
}
sub _createReadsPlacedFile($$){
  my($self,$args)=@_;
  my $outputFilename="$$args{baseDir}/reads.placed";
  
  my $read454Reads=$self->_read454ReadStatusFile("$$args{newblerDir}/454ReadStatus.txt");
  my $trim454Reads=$self->_read454TrimStatusFile("$$args{newblerDir}/454TrimStatus.txt");

  # merge hashes into readHash
  my($readHash);
  for my $name (keys (%$read454Reads)){
    while(my($k,$v)=each(%{ $$read454Reads{$name} })){ $$readHash{$name}{$k}=$v; }
    while(my($k,$v)=each(%{ $$trim454Reads{$name} })){ $$readHash{$name}{$k}=$v; }
  }

  # reads.placed - has the following fields
  #   NCBI TI number (just put * for unknown)
  #   read name
  #   start of trimmed read (1-based coordinates in original-read-space)
  #   number of bases in trimmed read
  #   orientation on contig (0=forward, 1=reverse)
  #   contig name
  #   supercontig name
  #   approximate start of trimmed read on contig (1-based coordinates)
  #   approximate start of trimmed read on supercontig (1-based coordinates)
  open(OUT,">",$outputFilename) or die "Could not open $outputFilename for writing";
  while (my ($name,$line)=each(%$readHash)){
    next if($$line{status}!~/assembl/i); # skip if not assembled or partially assembled
    my $superContigName='*';
    my $superContigStart='*';
    print OUT join("\t",($$line{TI},$$line{name},$$line{trimStart},$$line{trimmedLength},$$line{orientation},$$line{contig5Name},$superContigName,$$line{contig5Pos},$superContigStart));
    print OUT "\n";
  }
  close OUT;

  return $outputFilename;
}

sub _read454ReadStatusFile{
  my($self,$path)=@_;
  my $content;
  open(FILE,$path) or die "Could not open the file $path for reading";
  <FILE>; # burn the header
  while(<FILE>){
    my(%line);
    chomp;
    ($line{name},$line{status},$line{contig5Name},$line{contig5Pos},$line{contig5Strand},$line{contig3Name},$line{contig3Pos},$line{contig3Strand})=split /\t/;
    $line{orientation}='0';
    $line{orientation}='1' if($line{contig5Strand} eq '-');
    $line{TI}=$line{name};
    $$content{$line{name}}=\%line;
  }
  close FILE;
  return $content;
}
sub _read454TrimStatusFile{
  my($self,$path)=@_;
  my $content;
  open(FILE,$path) or die "Could not open the file $path for reading";
  <FILE>; # burn the header
  while(<FILE>){
    my(%line);
    chomp;
    ($line{name},$line{trimPoints},$line{trimmedLength},$line{origTrimpoints},$line{origLength},$line{rawLength})=split /\t/;
    ($line{trimStart},$line{trimStop})=split(/\-/,$line{trimPoints});
    $$content{$line{name}}=\%line;
  }
  close FILE;
  return $content;
}

1; # Oh why Perl, why must you have this true statement at the end?


