#!/usr/bin/env perl

# run_assembly_readMetrics.pl: puts out metrics for a raw reads file
# Author: Lee Katz <lkatz@cdc.gov>

package PipelineRunner;
my ($VERSION) = ('$Id: $' =~ /,v\s+(\d+\S+)/o);

my $settings = {
    appname => 'cgpipeline',
    # these are the subroutines for all read metrics
    metrics=>[qw(avgReadLength totalBases maxReadLength minReadLength avgQuality numReads)],
    #metrics=>[qw(avgReadLength totalBases maxReadLength minReadLength avgQuality avgQualPerRead numReads)],
    # these are the subroutines for all standard read metrics
    #stdMetrics=>[qw(N50 genomeLength numContigs assemblyScore)],
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

my @fastaExt=qw(.fasta .fa .mfa .fas .fna .sff.fna .sff.fasta);
my @fastqExt=qw(.fastq .fq .fastq.gz .fq.gz);
my @sffExt=qw(.sff);
$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH = $FSFind::LOG || *STDERR; print $FH "$0: ".(caller(1))[3].": @_\n";}

exit(main());

sub main() {
  $settings = AKUtils::loadConfig($settings);
  die(usage($settings)) if @ARGV<1;

  my @cmd_options=qw(help fast qual_offset=i minLength=i);
  GetOptions($settings, @cmd_options) or die;
  $$settings{qual_offset}||=33;
  $$settings{numcpus}||=AKUtils::getNumCPUs();
  $$settings{tempdir}||=AKUtils::mktempdir();

  # TODO read all files in however into a queue of reads.  Then streamline the analysis of those reads.
  my %final_metrics;
  print join("\t",qw(File avgReadLength totalBases maxReadLength minReadLength avgQuality numReads readScore))."\n";
  for my $input_file(@ARGV){
    my($filename,$dirname,$ext)=fileparse($input_file,(@sffExt, @fastaExt,@fastqExt));
    
    # fast settings: do a word count and multiply by the first few read lengths' average
    if($$settings{fast} && grep(/$ext/,@fastqExt)){
      $$settings{is_fastqGz}=1 if($ext=~/\.gz/);
      fastqFastStats($input_file,$settings);
    } elsif(grep(/$ext/,@fastqExt)){
      $$settings{is_fastqGz}=1 if($ext=~/\.gz/);
      my $entryQueue=Thread::Queue->new();    # for storing fastq 4-line entries
      my $metricsQueue=Thread::Queue->new(); # for receiving read lengths and metrics
      my $file=File::Spec->rel2abs($input_file);
      die("Input or file $file not found") unless -f $file;
      my @thr;
      $thr[$_]=threads->new(\&fastqIndividualReadMetrics,$entryQueue,$metricsQueue,$settings) for(0..$$settings{numcpus} - 1);
      readFastq($input_file,$entryQueue,$settings);
      $entryQueue->enqueue(undef) for(@thr);
      $_->join for(@thr);
      $metricsQueue->enqueue(undef); # this kills the crab
      fastqStats($input_file,$metricsQueue,$settings);
    } elsif(grep(/$ext$/,@fastaExt)){
      fastaStats($input_file,$settings);
    } elsif(grep(/$ext/,@sffExt)){
      sffStats($input_file,$settings);
    } else {
      logmsg "WARNING: I do not understand the extension $ext in the filename $filename.  SKIPPING.";
    }
  }

  return 0;
}

sub fastqFastStats{
  my($infile,$settings)=@_;
  my $qual_offset=$$settings{qual_offset};
  my $numReadsToAverage||=5000;
  my $numLinesToAverage=$numReadsToAverage*4; # calc it once here to save on cpu time
  my $firstReadlengthsTotalSize=0;
  my $firstReadlengthsTotalQual=0;
  my $numLines=0;
  if($$settings{is_fastqGz}){
    open(FASTQ,"gunzip -c $infile|") or die "Could not open the fastq $infile because $!";
  } else {
    open(FASTQ,"<",$infile) or die "Could not open the fastq $infile because $!";
  }

  # count all lines and pull metrics from the first XXXX reads
  while(<FASTQ>){
    chomp;
    $numLines++;
    next if($numLinesToAverage < $numLines);

    # 1=> readid, 2=>sequence,...
    my $mod=$numLines % 4;
    if($mod==2){
      $firstReadlengthsTotalSize+=length($_);
    }elsif($mod==0){
      my $qual=$_;
      my @qual=map(ord($_)-$qual_offset, split(//,$qual));
      $firstReadlengthsTotalQual+=sum(@qual);
    }
  }
  close FASTQ;

  my $avgReadLength=$firstReadlengthsTotalSize/$numReadsToAverage;
  my $avgQuality=$firstReadlengthsTotalQual/$firstReadlengthsTotalSize;
  my $numReads=$numLines/4;
  my $totalBases=$avgReadLength*$numReads;
  my $readScore=round(log($totalBases*$avgQuality*$avgReadLength));

  print join("\t",$infile,round($avgReadLength),int($totalBases),'.','.',round($avgQuality),$numReads,$readScore)."\n";
  return 1;
}

sub fastaStats{
  my($infile,$settings)=@_;
  my ($filename,$dirname,$ext)=fileparse($infile,@fastaExt);
  my $basename="$dirname/$filename";
  my $seqs=AKUtils::readMfa($infile,$settings);
  my $maxReadLength=1;
  my $minReadLength=99999999999999;
  my $readScore='.';
  my $totalBases=0;
  my $numReads=scalar(values(%$seqs));
  for(values(%$seqs)){
    my $length=length($_);
    $maxReadLength=$length if($length>$maxReadLength);
    $minReadLength=$length if($length<$minReadLength);
    $totalBases+=$length;
  }
  my $avgReadLength=round($totalBases/$numReads);

  # can I find a qual file to give quality scores?
  my $qualfile;
  $qualfile="$infile.qual" if(-e "$infile.qual");
  $qualfile="$basename.qual" if(-e "$basename.qual");
  my $avgQuality=".";
  my $readScore=round(log($totalBases*$avgReadLength*40)); # 40 gives a kind of stand-in quality score if there is none
  if($qualfile){
    my $totalQuality=0;
    my $quals=AKUtils::readMfa($qualfile,{keep_whitespace=>1});
    while(my($seqid,$qual)=each(%$quals)){
      $totalQuality+=$_ for(split(/\s+/,$qual));
    }
    $avgQuality=round($totalQuality/$totalBases);
    $readScore=round(log($totalBases*$avgQuality*$avgReadLength));
  }

  print join("\t",$infile,$avgReadLength,$totalBases,$maxReadLength,$minReadLength,$avgQuality,$numReads,$readScore)."\n";
  return 1;
}

sub sffStats{
  my($infile,$settings)=@_;
  system("sffinfo -s $infile > '$$settings{tempdir}/$infile.fna'");
  system("sffinfo -q $infile > '$$settings{tempdir}/$infile.fna.qual'");
  return fastaStats("$$settings{tempdir}/$infile.fna",$settings);
}

sub fastqStats{
  my($infile,$metricsQueue,$settings)=@_;
  my(@length,@qual);
  my $maxReadLength=1;
  my $minReadLength=99999999999999;
  my $readScore='.';
  while(defined(my $stats=$metricsQueue->dequeue)){
    my($length,$qual)=split(/_/,$stats);
    push(@length,$length);
    push(@qual,$qual);

    $minReadLength=$length if($minReadLength>$length);
    $maxReadLength=$length if($maxReadLength<$length);
  }
  my $totalBases=sum(@length);
  my $numReads=@length;
  my($avgQuality,$avgReadLength)=qw(. .);
  if($totalBases>1){
    my $totalQuality=0;
    for(my $i=0;$i<$numReads;$i++){
      $totalQuality+=$qual[$i]*$length[$i];
    }
    $avgQuality=round($totalQuality/$totalBases);
    $avgReadLength=round($totalBases/$numReads);
    $readScore=round(log($totalBases*$avgQuality*$avgReadLength));
    #my $readsScore=log($$m{totalBases} * $avgQuality * $avgReadLength);
  }

  print join("\t",$infile,$avgReadLength, $totalBases, $maxReadLength, $minReadLength, $avgQuality, $numReads, $readScore)."\n";
  
  return 1;
}

sub fastqIndividualReadMetrics{
  my($entryQueue,$metricsQueue,$settings)=@_;
  my $qual_offset=$$settings{qual_offset};
  my $i=0;
  my @cache;
  while(defined(my $fastqEntry=$entryQueue->dequeue)){
    my($id,$sequence,$plus,$qual)=split(/\n/,$fastqEntry);
    # get read length and quality
    my $length=length($sequence);
    my @qual=map(ord($_)-$qual_offset, split(//,$qual));
    my $avgQuality=sum(@qual)/@qual;
    push(@cache,join("_",$length,$avgQuality));
    if(@cache>1000000){
      $metricsQueue->enqueue(@cache);
      @cache=();
    }
  }
  $metricsQueue->enqueue(@cache);
  return $i;
}


sub readFastq{
  my($fastq,$entryQueue,$settings)=@_;
  my $i=0;
  if($$settings{is_fastqGz}){
    open(FASTQ,"gunzip -c $fastq|") or die "Could not open the fastq $fastq because $!";
  } else {
    open(FASTQ,"<",$fastq) or die "Could not open the fastq $fastq because $!";
  }
  my @cache;
  while(my $entry=<FASTQ>){
    $entry.=<FASTQ> for(1..3);
    push(@cache,$entry);
    if((++$i % ($$settings{numcpus}*100000)) == 0){
      $entryQueue->enqueue(@cache);
      @cache=();
      while($entryQueue->pending>10000000){ # 10 mill
        print $entryQueue->pending.".";
        sleep 1;
      }
    }
  }
  $entryQueue->enqueue(@cache);
  close FASTQ;
  return 1;
}

sub printMetrics{
  my ($metrics,$settings)=@_;
  my $header=$$settings{metrics};
  my $d="\t"; # field delimiter
  push(@$header,"readsScore");

  print "File$d".join($d,@$header)."\n";
  while(my($file,$m)=each(%$metrics)){
    my $avgQuality=$$m{avgQuality} || 40;
    my $avgReadLength=$$m{avgReadLength} || 100;
    my $readsScore=log($$m{totalBases} * $avgQuality * $avgReadLength);
    $$m{readsScore}=$readsScore;

    my @line;
    for(@$header){
      push(@line,$$m{$_} || ".");
    }

    print $file.$d.join($d,@line)."\n";
  }
  return 1;
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
    A reads file can be fasta or fastq
    The quality file for a fasta file is assumed to be reads.fasta.qual
  --fast for fast mode: fewer stats but much faster!
  --qual_offset 33
    Set the quality score offset (usually it's 33, so the default is 33)
  --minLength 0
    Set the minimum read length used for calculations
  --help for this help menu
  The reads score is log(numBases*avgQuality*avgReadLength), and 40 is the default quality and 100 is the default read length, if not shown
  "
}
