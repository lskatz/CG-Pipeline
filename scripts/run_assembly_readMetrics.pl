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
use Statistics::Descriptive;

use threads;
use Thread::Queue;

my @fastaExt=qw(.fasta .fa .mfa .fas .fna);
my @fastqExt=qw(.fastq .fq .fastq.gz .fq.gz);
my @sffExt=qw(.sff);
my @samExt=qw(.sam .bam);
$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH = $FSFind::LOG || *STDERR; print $FH "$0: ".(caller(1))[3].": @_\n";}

exit(main());

sub main() {
  $settings = AKUtils::loadConfig($settings);
  die(usage($settings)) if @ARGV<1;

  my @cmd_options=qw(help fast qual_offset=i minLength=i numcpus=i expectedGenomeSize=s histogram|hist reads-per-qual tempdir=s);
  GetOptions($settings, @cmd_options) or die $!;
  die usage() if($$settings{help});
  die "ERROR: need reads file\n".usage() if(@ARGV<1);
  $$settings{qual_offset}||=33;
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=AKUtils::mktempdir();
  $$settings{bufferSize}||=100000;
  $$settings{minLength}||=1;
  # the sample frequency is 100% by default or 1% if "fast"
  $$settings{sampleFrequency} ||= ($$settings{fast})?0.01:1;

  # get the right headers
  if($$settings{histogram}){
    print join("\t",qw(ReadLength Count Freq))."\n";
  } elsif($$settings{'reads-per-qual'}){
    print join("\t",qw(Phred Count))."\n";
  } else {
    print join("\t",qw(File avgReadLength totalBases minReadLength maxReadLength avgQuality numReads PE? coverage readScore medianFragmentLength))."\n";
  }

  for my $input_file(@ARGV){
    printReadMetricsFromFile($input_file,$settings);
  }

  return 0;
}

# main subroutine to print the metrics from a raw reads file
sub printReadMetricsFromFile{
  my($file,$settings)=@_;
  my($basename,$dirname,$ext)=fileparse($file,(@fastaExt,@fastqExt, @sffExt, @samExt));
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
  } elsif(grep(/$ext/,@samExt)) {
    $numEntries=readSam($file,$Q,$settings);
  } else {
    die "Could not understand filetype $ext";
  }

  # Avoid zero-read errors
  if($numEntries<1){
    logmsg "WARNING: there were no reads in $file. Moving on...\n";
    next;
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
    push(@{$count{tlen}},@{$$c{tlen}});
    push(@{$count{readLength}},@{$$c{readLength}});
    push(@{$count{readQual}},@{$$c{readQual}});
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
  my $medianFragLen='.';
  if(grep(/$ext/,@samExt)){
    # Bam files are PE if they have at least some fragment sizes.
    # See if tlen is present so that it can be calculated.
    $isPE=(@{$count{tlen}} > 10)?"yes":"no"; 
    if($isPE eq 'yes'){
      my $tlenStats=Statistics::Descriptive::Full->new;
      $tlenStats->add_data(@{$count{tlen}});
      my $tlen25=$tlenStats->percentile(25);
      my $tlen50=$tlenStats->percentile(50);
      my $tlen75=$tlenStats->percentile(75);
      $medianFragLen="$tlen50\[$tlen25,$tlen75]";
    }
  }
  # coverage is bases divided by the genome size
  my $coverage=($$settings{expectedGenomeSize})?round($count{extrapolatedNumBases}/$$settings{expectedGenomeSize}):'.';

  if($$settings{histogram}){
    printHistogram($count{readLength},$fractionReadsRead,$settings);
  } elsif($$settings{'reads-per-qual'}){
    printHistogram($count{readQual},$fractionReadsRead,$settings);
  } else {
    print join("\t",$file,$avgReadLength,$count{extrapolatedNumBases},$count{minReadLength},$count{maxReadLength},$avgQual,$count{numReads},$isPE,$coverage,'.',$medianFragLen)."\n";
  }
  return \%count;
}

sub printHistogram{
  my($data,$coefficient,$settings)=@_;
  my $numData=@$data;
  my %count;
  for(@$data){
    $_=sprintf("%0.0f",$_);
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
  my @readQual;
  my @tlen; # fragment length
  while(defined(my $tmp=$Q->dequeue)){
    my($seq,$qual,$tlen)=@$tmp;
    # trim and chomp
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

    $tlen||=0;
    push(@tlen,$tlen) if($tlen>0);

    # quality metrics
    my @qual;
    if($qual=~/\s/){ # if it is numbers separated by spaces
      @qual=split /\s+/,$qual;
    } else {         # otherwise, encoded quality
      @qual=map(ord($_)-$qual_offset, split(//,$qual));
    }
    $count{qualSum}+=sum(@qual);
    #push(@readQual, sum(@qual)/scalar(@qual));
    push(@readQual, sum(@qual));
  }
  $count{minReadLength}=$minReadLength;
  $count{maxReadLength}=$maxReadLength;
  $count{readLength}=\@length;
  $count{tlen}=\@tlen;
  $count{readQual}=\@readQual;
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

  # read the first one so that there is definitely going to be at least one read in the results
  my $bufferSize=$$settings{bufferSize};
  <$fp>; my $firstSeq=<$fp>; <$fp>; my $firstQual=<$fp>;
  my @queueBuffer=([$firstSeq,$firstQual]);
  my $numEntries=1;
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
  my $bufferSize=$$settings{bufferSize};
  <FNA>; my $firstSeq=<FNA>; <QUAL>; my $firstQual=<QUAL>;
  my @queueBuffer=([$firstSeq,$firstQual]);
  my $numEntries=1;
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

sub readSam{
  my($file,$Q,$settings)=@_;
  my($basename,$dir,$ext)=fileparse($file,@samExt);
  if($ext=~/sam/){
    open(SAM,$file) or die "ERROR: I could not read $file: $!";
  } elsif($ext=~/bam/){
    open(SAM,"samtools view $file | ") or die "ERROR: I could not use samtools to read $file: $!";
  } else {
    die "ERROR: I do not know how to read the $ext extension in $file";
  }

  my $bufferSize=$$settings{bufferSize};
  my $firstLine=<SAM>;
  my(undef,undef,undef,undef,undef,undef,undef,undef,$firstTlen,$firstSeq,$firstQual)=split /\t/,$firstLine;
  my @queueBuffer=([$firstSeq,$firstQual,$firstTlen]);
  my $numEntries=1;
  while(<SAM>){
    next if(/^@/);
    chomp;
    my($qname,$flag,$rname,$pos,$mapq,$cigar,$rnext,$pnext,$tlen,$seq,$qual)=split /\t/;
    #push(@queueBuffer,[$seq,$qual,$tlen]) if(rand() <= 0.1);
    push(@queueBuffer,[$seq,$qual,$tlen]) if(rand() <= $$settings{sampleFrequency});
    next if(++$numEntries % $bufferSize !=0);
    #print Dumper \@queueBuffer;die;
    
    $Q->enqueue(@queueBuffer);
    @queueBuffer=();

    while($Q->pending > $bufferSize * 3){
      sleep 1;
    }
  }
  $Q->enqueue(@queueBuffer);
  close SAM;
  
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
  "Prints useful read statistics on one or more read sets
  Usage: $0 *.fastq.gz
         $0 reads.fastq.gz | column -t
    A reads file can be fasta, sff, fastq, or fastq.gz
    The quality file for a fasta file reads.fasta is assumed to be reads.fasta.qual
  
  --fast            fast mode: samples 1% of the reads and extrapolates
  --numcpus     1   Number of cpus
  --qual_offset 33  Set the quality score offset
  --minLength   1   Set the minimum read length used for calculations
  -e        4000000 expected genome size, in bp
  
  HISTOGRAM
  --hist            Generate a histogram of the reads per length
  --reads-per-qual
  "
}
