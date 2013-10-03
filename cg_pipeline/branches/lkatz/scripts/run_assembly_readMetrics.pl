#!/usr/bin/env perl

# run_assembly_readMetrics.pl: puts out metrics for a raw reads file
# Author: Lee Katz <lkatz@cdc.gov>

package PipelineRunner;
my ($VERSION) = ('$Id: $' =~ /,v\s+(\d+\S+)/o);

my $settings = {
    appname => 'cgpipeline',
    # these are the subroutines for all read metrics
    metrics=>[qw(avgReadLength totalBases maxReadLength minReadLength avgQuality numReads)],
};

use strict;
no strict "refs";
use FindBin;
use lib "$FindBin::RealBin/../lib";
$ENV{PATH} = "$FindBin::RealBin:".$ENV{PATH};
use AKUtils qw(logmsg);

use Getopt::Long;
use File::Temp ('tempdir');
use File::Path;
use File::Spec;
use File::Copy;
use File::Basename;
use List::Util qw(min max sum shuffle);
use CGPipelineUtils;
use Data::Dumper;

use threads;
use Thread::Queue;

my @fastaExt=qw(.fasta .fa .mfa .fas .fna);
my @fastqExt=qw(.fastq .fq .fastq.gz .fq.gz);
my @sffExt=qw(.sff);
$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH = $FSFind::LOG || *STDERR; print $FH "$0: ".(caller(1))[3].": @_\n";}

exit(main());

sub main() {
  $settings = AKUtils::loadConfig($settings);
  die(usage($settings)) if @ARGV<1;

  my @cmd_options=qw(help fast qual_offset=i minLength=i numcpus=i expectedGenomeSize=s histogram);
  GetOptions($settings, @cmd_options) or die;
  die usage() if($$settings{help});
  die "ERROR: need reads file\n".usage() if(@ARGV<1);
  $$settings{qual_offset}||=33;
  $$settings{numcpus}||=AKUtils::getNumCPUs();
  $$settings{tempdir}||=AKUtils::mktempdir();
  $$settings{bufferSize}||=100000;
  $$settings{minLength}||=1;
  # the sample frequency is 100% by default or 1% if "fast"
  $$settings{sampleFrequency} ||= ($$settings{fast})?0.01:1;

  print join("\t",qw(File avgReadLength totalBases minReadLength maxReadLength avgQuality numReads PE? coverage readScore))."\n";
  for my $input_file(@ARGV){
    printReadMetricsFromFile($input_file,$settings);
  }

  return 0;
}

# main subroutine to print the metrics from a raw reads file
sub printReadMetricsFromFile{
  my($file,$settings)=@_;
  my($basename,$dirname,$ext)=fileparse($file,(@fastaExt,@fastqExt, @sffExt));
  # start the queue and threads
  my $Q=Thread::Queue->new();
  my @thr;
  $thr[$_]=threads->new(\&readMetrics,$Q,$settings) for(0..$$settings{numcpus}-1);

  my $numEntries=0;
  if(grep(/$ext/,@fastqExt)){
    $numEntries=readFastq($file,$Q,$settings);
  } elsif(grep(/$ext/,@fastaExt)) {
    $numEntries=readFasta($file,$Q,$settings);
  } elsif(grep(/$ext/,@sffExt)){
    $numEntries=readSff($file,$Q,$settings);
  } else {
    die "Could not understand filetype $ext";
  }

  # Combine the threads
  my %count=(minReadLength=>1e999); # need a min length to avoid a bug later
  $Q->enqueue(undef) for(@thr);
  for(@thr){
    my $c=$_->join;
    $count{numBases}+=$$c{numBases};
    $count{numReads}+=$$c{numReads};
    $count{qualSum} +=$$c{qualSum};
    $count{maxReadLength}=max($$c{maxReadLength},$count{maxReadLength});
    $count{minReadLength}=min($$c{minReadLength},$count{minReadLength});
    push(@{$count{readLength}},@{$$c{readLength}});
  }

  # extrapolate the counts to the total number of reads if --fast
  my $fractionReadsRead=$count{numReads}/$numEntries;
  $count{numReads}=$numEntries;
  $count{extrapolatedNumBases}=int($count{numBases}/$fractionReadsRead);
  $count{extrapolatedNumReads}=int($count{numReads}*$fractionReadsRead);

  # derive some more values
  my $avgQual=round($count{qualSum}/$count{numBases});
  my $avgReadLength=round($count{numBases}/$count{extrapolatedNumReads});
  my $isPE=(AKUtils::is_fastqPE($file))?"yes":"no";
  # coverage is bases divided by the genome size
  my $coverage=($$settings{expectedGenomeSize})?round($count{extrapolatedNumBases}/$$settings{expectedGenomeSize}):'.';

  print join("\t",$file,$avgReadLength,$count{extrapolatedNumBases},$count{minReadLength},$count{maxReadLength},$avgQual,$count{numReads},$isPE,$coverage,'.')."\n";
  printHistogram($count{readLength},$fractionReadsRead,$settings) if($$settings{histogram});
  return \%count;
}

sub printHistogram{
  my($data,$coefficient,$settings)=@_;
  my @data=map(int($_/100)*100,@$data);
  my $numData=@data;
  my %count;
  for(@data){
    $count{$_}++;
  }
  #$sum=int($sum / $coefficient);
  #print Dumper \%count;die;
  for my $datum(sort {$a<=>$b} keys(%count)){
    print join("\t",$datum,$count{$datum},round($count{$datum}/$numData))."\n";
  }
  print join("\t","total",$numData,'.')."\n";
  return \%count;
}

# Reads a Thread::Queue to give metrics but does not derive any metrics, e.g. avg quality
sub readMetrics{
  my($Q,$settings)=@_;
  my $qual_offset=$$settings{qual_offset} || die "Internal error";
  my %count;
  my $minReadLength=1e999;
  my $maxReadLength=0;
  my @length;
  while(defined(my $tmp=$Q->dequeue)){
    my($seq,$qual)=@$tmp;
    # trim
    $seq =~s/^\s+|\s+$//g;
    $qual=~s/^\s+|\s+$//g;
    my $readLength=length($seq);
    next if($readLength<$$settings{minLength});
    push(@length,$readLength);
    $count{numBases}+=$readLength;
    if($readLength<$minReadLength){
      $minReadLength=$readLength;
    } elsif ($readLength>$maxReadLength){
      $maxReadLength=$readLength;
    }
    $count{numReads}++;

    # quality metrics
    my @qual;
    if($qual=~/\s/){ # if it is numbers separated by spaces
      @qual=split /\s+/,$qual;
    } else {         # otherwise, encoded quality
      @qual=map(ord($_)-$qual_offset, split(//,$qual));
    }
    $count{qualSum}+=sum(@qual);
  }
  $count{minReadLength}=$minReadLength;
  $count{maxReadLength}=$maxReadLength;
  $count{readLength}=\@length;
  return \%count;
}

sub readFastq{
  my($file,$Q,$settings)=@_;
  my($basename,$dirname,$ext)=fileparse($file,(@fastaExt,@fastqExt, @sffExt));
  my $fp;
  if($ext=~/\.gz/){
    open($fp,"gunzip -c $file |") or die "Could not open $file:$!";
  } else {
    open($fp,$file) or die "Could not open fastq $file:$!";
  }
  my @queueBuffer;
  my $numEntries=0;
  my $bufferSize=$$settings{bufferSize};
  while(<$fp>){
    $numEntries++;
    my $seq=<$fp>;
    <$fp>; # burn the "plus" line
    my $qual=<$fp>;
    push(@queueBuffer,[$seq,$qual]) if(rand() <= $$settings{sampleFrequency});
    next if($numEntries % $bufferSize !=0);
    # Don't process the buffer until it is full

    # flush the buffer
    $Q->enqueue(@queueBuffer);
    @queueBuffer=();
    # pause if the queue is too full
    while($Q->pending > $bufferSize * 3){
      sleep 1;
    }
    if($Q->pending < $bufferSize && $numEntries > $bufferSize){
      #logmsg "The queue is getting emptied too fast. Increase the buffer size to optimize this script. Queue is at ".$Q->pending;
      $bufferSize*=2; # double the buffer size then
      #logmsg "Buffer size is now at $bufferSize";
    }
  }
  $Q->enqueue(@queueBuffer);
  close $fp;
  return $numEntries;
}

sub readSff{
  my($file,$Q,$settings)=@_;
  my($basename,$dirname,$ext)=fileparse($file,(@fastaExt,@fastqExt, @sffExt));

  local $/="\n>";
  open(FNA,"sffinfo -s $file | ") or die "Could not open $file:$!";
  open(QUAL,"sffinfo -q $file | ") or die "Could not open $file:$!";
  my @queueBuffer;
  my $numEntries=0;
  my $bufferSize=$$settings{bufferSize};
  while(my $defline=<FNA>){
    $numEntries++;
    <QUAL>; # burn the qual defline because it is the same as the fna
    my $seq=<FNA>;
    my $qual=<QUAL>;
    push(@queueBuffer,[$seq,$qual]);
    next if($numEntries % $bufferSize !=0);
    # Don't process the buffer until it is full

    # flush the buffer
    $Q->enqueue(@queueBuffer);
    @queueBuffer=();
    if($$settings{fast} && $numEntries>$bufferSize){
      while(<FNA>){
        $numEntries++; # count the rest of the reads
      }
      last;
    }
    # pause if the queue is too full
    while($Q->pending > $bufferSize * 3){
      sleep 1;
    }
  }
  $Q->enqueue(@queueBuffer);
  close QUAL; close FNA;

  return $numEntries;
}

sub readFasta{
  my($file,$Q,$settings)=@_;
  my($basename,$dirname,$ext)=fileparse($file,(@fastaExt,@fastqExt, @sffExt));

  local $/="\n>";
  open(FNA,$file) or die "Could not open $file:$!";
  open(QUAL,"$file.qual") or warn "WARNING: Could not open qual $file.qual:$!";
  my @queueBuffer;
  my $numEntries=0;
  my $bufferSize=$$settings{bufferSize};
  while(my $defline=<FNA>){
    $numEntries++;
    <QUAL>; # burn the qual defline because it is the same as the fna
    my $seq=<FNA>;
    my $qual=<QUAL> || "";
    push(@queueBuffer,[$seq,$qual]);
    next if($numEntries % $bufferSize !=0);
    # Don't process the buffer until it is full

    # flush the buffer
    $Q->enqueue(@queueBuffer);
    @queueBuffer=();
    if($$settings{fast} && $numEntries>$bufferSize){
      while(<FNA>){
        $numEntries++; # count the rest of the reads
      }
      last;
    }
    # pause if the queue is too full
    while($Q->pending > $bufferSize * 3){
      sleep 1;
    }
  }
  $Q->enqueue(@queueBuffer);
  close FNA; close QUAL;

  return $numEntries;
}

# Truncate to the hundreds place.
# Yes I understand it's not technically rounding.
sub round{
  my ($num)=(@_);
  my $rounded=int($num*100)/100;
  return sprintf("%.2f",$rounded); # put in zero padding in case it truncates at a zero in the hundreds place
}

sub usage{
  my ($settings)=@_;
  "Prints useful assembly statistics
  Usage: $0 reads.fasta 
         $0 reads.fasta | column -t
    A reads file can be fasta, sff, or fastq
    The quality file for a fasta file reads.fasta is assumed to be reads.fasta.qual
  --fast for fast mode: samples 1% of the reads and extrapolates
  -n 1 to specify the number of cpus (default: all cpus)
  --qual_offset 33
    Set the quality score offset (usually it's 33, so the default is 33)
  --minLength 1
    Set the minimum read length used for calculations
  -e 4000000 expected genome size, in bp
  --hist to generate a histogram of the reads
  "
}
