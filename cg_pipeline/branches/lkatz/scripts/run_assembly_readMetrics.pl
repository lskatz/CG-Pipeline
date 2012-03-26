#!/usr/bin/env perl

# run_assembly_readMetrics.pl: puts out metrics for a raw reads file
# Author: Lee Katz <lkatz@cdc.gov>

package PipelineRunner;
my ($VERSION) = ('$Id: $' =~ /,v\s+(\d+\S+)/o);

my $settings = {
    appname => 'cgpipeline',
};
my $stats;

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

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH = $FSFind::LOG || *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

exit(main());

sub main() {
  $settings = AKUtils::loadConfig($settings);
  die(usage($settings)) if @ARGV<1;

  my @cmd_options=qw(help fast);
  GetOptions($settings, @cmd_options) or die;

  for my $input_file(@ARGV){
    my $file=File::Spec->rel2abs($input_file);
    die("Input or file $file not found") unless -f $file;
  
    my %metrics;
    if($$settings{fast}){
      %metrics=fastFastqStats($file,$settings);
    } else {
      %metrics=readMetrics($file,$settings);
    }
    print Dumper \%metrics;
  }
  return 0;
}

sub readMetrics{
  my($file,$settings)=@_;

  # STEP 1: READ THE FILE
  my($seqs,$qual);
  my $ext=(split(/\./,$file))[-1];
  my($name,$path,$ext)=fileparse($file,qw(.fastq.gz .fastq .fq .fa .fa .fas .fasta .fna .mfa));
  $$settings{is_fastqGz}=1 if($ext=~/fastq.gz/);
  if($ext=~/fastq|fq|fastq.gz/){
    ($seqs,$qual)=readFastq($file,$settings);
  }
  elsif($ext=~/fa|fas|fasta|fna|mfa/i){
    $seqs=AKUtils::readMfa($file,$settings);
    my $qualfile;
    for("$name$ext.qual","$name.qual"){
      $qualfile=$_ if(-e $_);
    }
    # convert the quality file to Illumina qual strings
    $qual=AKUtils::readMfa($qualfile,{keep_whitespace=>1});
    while(my($id,$qualstr)=each(%$qual)){
      my $newQual;
      my @tmpQual=split(/\s+/,$qualstr);
      $newQual.=chr($_+33) for(@tmpQual);
      $$qual{$id}=$newQual;
    }
  } else {
    die "$ext extension not supported";
  }
  
  # STEP 1b: FILL IN QUALITY VALUES FOR READS WHOSE QUALITY IS NULL
  if(!keys(%$qual)){
    warn "Warning: There are no quality scores associated with $file";
    my $nullQual=chr(33);
    $$qual{$_}=$nullQual x length($$seqs{$_}) for(keys(%$seqs));
  }

  # STEP 2: METRICS
  my $seqCounter=0;
  my($totalReadLength,$maxReadLength,$totalReadQuality,$totalQualScores,$avgReadQualTotal)=(0,0);
  while(my($id,$seq)=each(%$seqs)){
    # read metrics
    my $readLength=length($seq);
    $totalReadLength+=$readLength;
    $maxReadLength=$readLength if($readLength>$maxReadLength);
    
    # quality metrics
    my $qualStr=$$qual{$id};
    my @qual;
    if(!$qualStr){
      warn "Warning: Could not find qual for $id. Setting to 0.";
      $qualStr=chr(33)x$readLength;
    }
    @qual=map(ord($_)-33,split(//,$qualStr));
    my $sumReadQuality=sum(@qual);
    my $thisReadAvgQual=$sumReadQuality/$readLength;
    $avgReadQualTotal+=$thisReadAvgQual;
    $totalReadQuality+=$sumReadQuality;
    $totalQualScores+=$readLength;
    
    $seqCounter++;
  }
  my $avgReadLength=$totalReadLength/$seqCounter;
  my $avgQuality=$totalReadQuality/$totalQualScores;
  my $avgQualPerRead=$avgReadQualTotal/$seqCounter;

  my %metrics=(
    avgReadLength=>$avgReadLength,
    totalBases=>$totalReadLength,
    maxReadLength=>$maxReadLength,
    avgQuality=>$avgQuality,
    avgQualPerRead=>$avgQualPerRead,
  );
  return %metrics;
}

sub readFastq{
  my($fastq,$settings)=@_;
  my $seqs={};
  my $quals={};
  my $i=0;
  if($$settings{is_fastqGz}){
    open(FASTQ,"gunzip -c $fastq|") or die "Could not open the fastq $fastq because $!";
  } else {
    open(FASTQ,"<",$fastq) or die "Could not open the fastq $fastq because $!";
  }
  while(my $id=<FASTQ>){
    my $sequence=<FASTQ>;
    my $plus=<FASTQ>;
    my $qual=<FASTQ>;
    next if(!$sequence);
    $id=~s/^@//;
    $$seqs{$id}=$sequence;
    $$quals{$id}=$qual;
  }
  close FASTQ;
  return ($seqs,$quals);
}
sub fastFastqStats{
  my($fastq,$settings)=@_;
  my($name,$path,$ext)=fileparse($fastq,qw(.fastq.gz .fastq .fq .fa .fa .fas .fasta .fna .mfa));
  $$settings{is_fastqGz}=1 if($ext=~/fastq.gz/);
  my %metrics;
  my $command="wc -l $fastq";
    $command="gunzip -c $fastq | wc -l" if($$settings{is_fastqGz});
  my $lines=`$command`; die "problem with wc -l" if $?;
  $metrics{numReads}=$lines/4;
  
  # bases in the first read
  $command="head -2 $fastq|tail -1";
  $command="gunzip -c $fastq | head -2 |tail -1" if($$settings{is_fastqGz});
  my $firstLine=`$command`; die "Problem with head or tail" if $?;
  chomp($firstLine);
  my $basesInFirstLine=length($firstLine);
  $metrics{totalBases}=$metrics{numReads}*$basesInFirstLine;

  return %metrics;
}

sub usage{
  my ($settings)=@_;
  "Prints useful assembly statistics
  Usage: $0 reads.fasta
    A reads file can be fasta or fastq
    The quality file for a fasta file is assumed to be reads.fasta.qual
  --fast for fast mode: fewer stats but much faster!
  --help for this help menu
  "
}
