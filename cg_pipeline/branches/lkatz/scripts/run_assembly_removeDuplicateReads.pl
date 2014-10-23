#!/usr/bin/env perl 

# run_assembly_trimClean: trim and clean a set of raw reads
# Author: Lee Katz <lkatz@cdc.gov>

package PipelineRunner;
my ($VERSION) = ('$Id: $' =~ /,v\s+(\d+\S+)/o);

use strict;
use FindBin;
use lib "$FindBin::RealBin/../lib";
$ENV{PATH} = "$FindBin::RealBin:".$ENV{PATH};
use AKUtils;

use Getopt::Long;
use File::Temp ('tempdir');
use File::Path;
use File::Spec;
use File::Copy;
use File::Basename;
use List::Util qw(min max sum shuffle);
use POSIX;

use Data::Dumper;
use threads;
use Thread::Queue;
use threads::shared;

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {print STDERR "$0: ".(caller(1))[3].": @_\n";}

my %threadStatus:shared;
my @IOextensions=qw(.fastq .fastq.gz .fq);

exit(main());

sub main{
  my $settings = {
      appname => 'cgpipeline',
  };

  # get settings from cg_pipeline/conf/cgpipelinerc and ./cgpipelinerc
  $settings=AKUtils::loadConfig($settings);
  # get CLI flags into $settings with Getopt::Long
  GetOptions($settings,qw(help tempdir=s downsample=s length=i sizeTo=s stdin stdin-compressed));
  die usage() if($$settings{help});
  # additional settings using ||= operator
  $$settings{poly}||=1; # SE by default
  $$settings{tempdir}||=AKUtils::mktempdir(); # any temp file you make should go under the tempdir
  $$settings{downsample}||=1;
  die "ERROR: downsample should be a number and less than 1. User requested $$settings{downsample}.\n".usage() if($$settings{downsample}=~/[A-Za-z]/ || $$settings{downsample}>1);

  my $infile=[@ARGV];
  
  my @readArr;
  if(@ARGV && !$$settings{stdin} && !$$settings{'stdin-compressed'}){
    for(@$infile){
      my $reads=readFastq($_,$settings);
      push(@readArr,@$reads);
    }
  } elsif($$settings{stdin} || $$settings{'stdin-compressed'}){
    logmsg "Reading from stdin";
    @readArr=@{ readFastqStdin($settings) };
  } else {
    die "ERROR: did not detect any file!\n".usage();
  }

  removeDuplicateReads(\@readArr,$settings);

  logmsg "Done";
  return 0;
}

sub readFastq{
  my($infile,$settings)=@_;
  my $linker="~" x 4; # tildes are probably not in the sequence, so it's safe for putting into the hash ID
  my $linkerFirstChar=substr($linker,0,1);

  my ($basename,$dir,$inSuffix) = fileparse($infile,@IOextensions);

  my $poly=AKUtils::is_fastqPE($infile,$settings)+1;
  $$settings{poly}=$poly if($poly==1); # if you see SE anywhere, then it is safer to assume SE for the entirety

  if($inSuffix=~/\.fastq$/){
    open(FQ,'<',$infile) or die "Could not open $infile for reading: $!";
  }
  elsif($inSuffix=~/\.fastq\.gz/){
    open(FQ,"gunzip -c $infile | ") or die "Could not open $infile for reading: $!";
  }
  else{
    die "Could not determine the file type for reading based on your extension $inSuffix";
  }
  logmsg "Reading $infile with poly=$poly";
  my $i=0;
  my @reads=();
  while(my $id=<FQ>){
    my $sequence=<FQ>; chomp($sequence);
    my $discard=<FQ>;
    my $qual=<FQ>;
    my $read="$id$sequence\n+\n$qual";
    for(my $j=1;$j<$poly;$j++){
      my $id2=<FQ>;
      my $sequence2=<FQ>; chomp($sequence2);
      my $discard2=<FQ>;
      my $qual2=<FQ>;
      $read.="$id2$sequence2\n+\n$qual2";
    }

    push(@reads,$read);
  }
  close FQ;
  
  return \@reads;
}

sub readFastqStdin{
  my($settings)=@_;
  my $linker="~" x 4; # tildes are probably not in the sequence, so it's safe for putting into the hash ID
  my $linkerFirstChar=substr($linker,0,1);

  my $poly=2; # TODO customize this
  $$settings{poly}=$poly;
  logmsg "Reading stdin with poly=$poly";
  my $i=0;
  my @reads=();
  if($$settings{'stdin-compressed'}){
    open(FQ,"gunzip -c |") or die "Could not gunzip from stdin: $!";
  } else {
    *FQ=*STDIN; # copy the filehandle
  }
  while(my $id=<FQ>){
    my $sequence=<FQ>; chomp($sequence);
    my $discard=<FQ>;
    my $qual=<FQ>;
    my $read="$id$sequence\n+\n$qual";
    for(my $j=1;$j<$poly;$j++){
      my $id2=<FQ>;
      my $sequence2=<FQ>; chomp($sequence2);
      my $discard2=<FQ>;
      my $qual2=<FQ>;
      $read.="$id2$sequence2\n+\n$qual2";
    }

    push(@reads,$read);
  }
  close FQ;
  
  return \@reads;
}

# Return the best read in an array of reads
# TODO maybe make a consensus quality score instead of just getting the best representative
sub bestReadOld{
  my($reads,$settings)=@_;


  my $numReads=scalar(@$reads);
  return $$reads[0] if($numReads<2);
  @$reads=sort {
    my(undef,undef,undef,$qual1a,undef,undef,undef,$qual2a)=split(/\n/,$a);
    my(undef,undef,undef,$qual1b,undef,undef,undef,$qual2b)=split(/\n/,$b);
    $qual1a||="";
    $qual1b||="";
    $qual2a||="";
    $qual2b||="";

    # calculate total quality. 
    # Avg doesn't matter because we only care about reads that are the same sequence with differing qualities.
    # Different but adjacent sequences will not matter.
    # Qual offset doesn't matter because all reads will have the same offset.
    my $sumA=sum(map(ord($_),split(//,"$qual1a$qual2a")));
    my $sumB=sum(map(ord($_),split(//,"$qual1b$qual2b")));
    return $sumB <=> $sumA;
  } @$reads;
  return $$reads[0];
}

sub bestRead{
  my($reads,$settings)=@_;

  my $numReads=scalar(@$reads);
  # the first read will work if there's only one read
  return $$reads[0] if($numReads<2);

  # assume certain characteristics of the combined read from the first read
  my($id1,$seq1,undef,undef,$id2,$seq2,undef,undef)=split(/\n/,$$reads[0]);
  my ($length1,$length2)=(length($seq1),length($seq2));
  
  # assume that since two reads agree, then their quality is top-notch
  my $qual1=chr(40+33) x $length1;
  my $qual2=chr(40+33) x $length2;

  # return the actual read
  #die "\n".join("",@$reads)."\n$id1\n$seq1\n+\n$qual1\n$id2\n$seq2\n+\n$qual2\n";
  return "$id1\n$seq1\n+\n$qual1\n$id2\n$seq2\n+\n$qual2\n" if($seq2);
  return "$id1\n$seq1\n+\n$qual1\n";
}

sub removeDuplicateReads{
  my($reads,$settings)=@_;
  my $linker="~" x 4; # tildes are probably not in the sequence, so it's safe for putting into the hash ID
  my $linkerFirstChar=substr($linker,0,1);
  my $l=$$settings{length}||0;
  my $poly=$$settings{poly};

  # change the downsample parameter if settings{sizeTo} is set
  if($$settings{sizeTo}){
    $$settings{downsample}=findDownsamplingFromSizeto($reads,$$settings{sizeTo},$settings);
  }

  # sort reads by quality first so that the best quality is retained.
  #sortReads($reads,$settings);
  #die Dumper [$$reads[0],$$reads[100]];

  # Compare sameness reads but with differing qualities: bin the reads
  # Downsample within this loop so that abundance can factor in.
  logmsg "Binning reads";
  my $i=0;
  my %binnedRead;
  my $numReads=@$reads; # or readPairs
  for (my $i=0;$i<$numReads;$i++){
    my($id1,$seq1,undef,$qual1,$id2,$seq2,undef,$qual2)=split(/\n/,$$reads[$i]);
    $seq2||="";
    my $hashId="$seq1$linker$seq2";
    if($l){
      $hashId=~s/^(.{$l,$l}).*$linker(.+)/$1$linker$2/; # accept only X nucleotides from the front
      $hashId=~s/($linker.{$l,$l}).*($|$linker)/$1$2/g if($poly>1);
    }
    # downsample here
    next if(rand() > $$settings{downsample});
    push(@{ $binnedRead{$hashId} }, $$reads[$i]);
    print STDERR "." if(++$i % 10000 == 0);
  }
  print STDERR "\n";
  undef($reads); # free up some space

  # Print one sequence read per duplicate set
  logmsg "Choosing best read per bin";
  $i=0;
  for my $readArr(values(%binnedRead)){
    print bestRead($readArr,$settings);
    print STDERR "." if(++$i % 100000 == 0);
  }
  print STDERR "\n";

  return scalar(keys(%binnedRead));
}

sub findDownsamplingFromSizeto{
  my($reads,$sizeTo,$settings)=@_;

  my $bp=0;
  my $numReads=@$reads;
  for(my $i=0;$i<$numReads;$i++){
    my(undef,$seq1,undef,undef,$seq2)=split(/\n/,$$reads[$i]);
    $seq2||="";
    $bp+=length($seq1)+length($seq2);
  }
  my $downsample=$sizeTo/$bp;
  $downsample=1 if($downsample>1);
  logmsg "Downsampling as ".sprintf("%0.4f",$downsample)." of the original with requested size of $sizeTo";
  return $downsample;
}


sub usage{
  "Removes duplicate reads from a raw read set
   Usage: $0 read.fastq[.gz] > read.fastq
     --downsample 0.1    # only keep 10% of the reads
     -size 1000000       # downsample to 1Mb. Internally, a new --downsample parameter is calculated and will be reported in stderr
     --length 100        # only consider up to 100bp when deciding if a read is a duplicate. Default: no limit
     --stdin             # read a file as stdin (uncompressed)
     --stdin-compressed  # read a file as stdin (compressed)
  EXAMPLES
    $0 read.fastq[.gz] | gzip -c > read.fastq.gz
    cat read.fastq.gz | $0 --stdin-compressed | gzip -c > noDupes.fastq.gz
  ";
}

