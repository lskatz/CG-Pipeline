#!/usr/bin/env perl 

# run_assembly_trimClean: trim and clean a set of raw reads
# Author: Lee Katz <lkatz@cdc.gov>
# TODO: make a low-memory binning method, combining only 
# a few reads at a time and perhaps using temporary
# intermediate files.

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
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; my $caller=(caller(1))[3]; $caller=~s/::main//; die("$0: $caller: ".$e); };
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
  GetOptions($settings,qw(help tempdir=s downsample=f length=i sizeTo|size=s stdin stdin-compressed highest=i bin!)) or die $!;
  # additional settings using ||= operator
  $$settings{bin}//=1;
  $$settings{poly}||=1; # SE by default
  $$settings{tempdir}||=AKUtils::mktempdir(); # any temp file you make should go under the tempdir
  $$settings{downsample}||=1;
  $$settings{highest}||=40+33; # I is the 9th letter... don't want to go past Z
  my $infile=[@ARGV];

  die usage($settings) if($$settings{help} || !@$infile);
  die "ERROR: downsample parameter $$settings{downsample} is not between 0 and 1.\n".usage($settings) if($$settings{downsample}>1 || $$settings{downsample} < 0);

  if($$settings{bin}){
    logmsg "User supplied --bin; Binning reads.";
  } else {
    logmsg "User supplied --nobin; not binning reads.";
  }

  
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
    die "ERROR: did not detect any file!\n".usage($settings);
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

# Return the only read or if there are duplicates,
# calculate a new quality score based on consensus.
# The read ID will be the same as the first read in the array.
sub bestRead{
  my($reads,$settings)=@_;

  my $numReads=scalar(@$reads);
  # the first read will work if there's only one read
  return $$reads[0] if($numReads<2);

  # assume certain characteristics of the combined read from the first read
  my($id1,$seq1,undef,$firstQual1,$id2,$seq2,undef,$firstQual2)=split(/\n/,$$reads[0]);
  my ($length1,$length2)=(length($seq1),length($seq2));
  # @p1, @p2 contain the prob(error) for each base along the read.
  # Later it will get multiplied against the other probs found in the other reads
  my @p1=map(10**(-1*(ord($_)-33)/10),split(//,$firstQual1));
  my @p2=map(10**(-1*(ord($_)-33)/10),split(//,$firstQual2));

  # Calculate a new error probability for each base in the duplicated read.
  # Start with the second read because the first read was used to initialize @p1,@p2
  for(my $i=1;$i<$numReads;$i++){
    my(undef,undef,undef,$q1,undef,undef,undef,$q2)=split(/\n/,$$reads[$i]);
    # F/R quality scores are put into probabilities for this read
    my @pForward=map(10**(-1*(ord($_)-33)/10),split(//,$q1));
    my @pReverse=map(10**(-1*(ord($_)-33)/10),split(//,$q2));

    # probability of error = p1 * p2 * ...
    for(my $j=0;$j<$length1;$j++){
      $p1[$j]*=$pForward[$j];
    }
    for(my $j=0;$j<$length2;$j++){
      $p2[$j]*=$pReverse[$j];
    }
    
  }
  my $highestScore=$$settings{highest};
  my $highestChr = chr($highestScore);
  # back to base-33 scores but not chr-ed yet
  for(@p1,@p2){$_=0.00000001 if($_<=0); }
  my @qual1=map(int(-10*log($_)/log(10)+33),@p1);
  my @qual2=map(int(-10*log($_)/log(10)+33),@p2);
  # don't let anything go past the highest score
  for(@qual1,@qual2){
    $_=$highestScore if($_ > $highestScore); # compare vs highest possible score
    $_=chr($_); # can now chr now that the comparison has been done
  }
  # Back to a string
  my $qual1=join("",@qual1);
  my $qual2=join("",@qual2);
  
  # return the actual read
  #die "\n".join("",@$reads)."\n$id1\n$seq1\n+\n$qual1\n$id2\n$seq2\n+\n$qual2\n";

  # Tack on how many reads were combined
  if($numReads>1){
    $id1.=" CGP:combinedReads=$numReads";
    $id2.=" CGP:combinedReads=$numReads";
  }
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
    # If we're not binning/deduplicating, then just have a unique ID
    $hashId="nobinning$i" if(!$$settings{bin});
    # bin the reads
    push(@{ $binnedRead{$hashId} }, $$reads[$i]);
    print STDERR "." if($i % 100000 == 0 && $i>0);
  }
  print STDERR "\n";
  logmsg "Found ".scalar(keys(%binnedRead))." unique read or read pairs out of ".scalar(@$reads);
  undef($reads); # free up some space

  # Print one sequence read per duplicate set
  logmsg "Combining identical reads per bin and recalculating quality scores";
  $i=0;
  while(my($hashId,$readArr)=each(%binnedRead)){
    # downsample here
    next if(rand() > $$settings{downsample});
    # print the best read out of the duplicates
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
    my(undef,$seq1,undef,undef,undef,$seq2)=split(/\n/,$$reads[$i]);
    $seq2||="";
    $bp+=length($seq1)+length($seq2);
  }
  my $downsample=$sizeTo/$bp;
  $downsample=1 if($downsample>1);
  logmsg "Downsampling as ".sprintf("%0.4f",$downsample)." of the original with requested size of $sizeTo";
  return $downsample;
}


sub usage{
  my($settings)=@_;
  
  "Removes duplicate reads from a raw read set.

   Usage: $0 read.fastq[.gz] > read.fastq
     --[no]bin           If --bin, Combine identical reads.
                         Any two quality scores are combined as a result of -log(p * p) 
                         where p is the probability of getting an error and is fraction 
                         representation of the phred score.
                         If --nobin, Don't deduplicate; You are using this
                         script purely for the downsampling.
     --downsample        Only keep a percentage of reads [default: 1.0]
     --size              Downsample to N base pairs. Internally, a new 
                         --downsample parameter is calculated and will 
                         be reported in stderr
     --length            Only consider up to X bp when deciding if a read 
                         is a duplicate. [Default: 0 (no limit)]
                         Warning: you might lose some sequence information 
                         when reads are binned if you use --length
     --stdin             Read a file as stdin (uncompressed)
     --stdin-compressed  Read a file as stdin (compressed)
     --highest           The highest allowed quality score after combining 
                         reads [Default: 73] The highest Illumina score is 
                         normally 73, which is 'I'. However if you want 'Z',
                         it is ".(40+33+(26-9)).".
  EXAMPLES
    $0 read.fastq[.gz] | gzip -c > read.fastq.gz
    $0 --stdin-compressed < read.fastq.gz | gzip -c > noDupes.fastq.gz
  ";
}

