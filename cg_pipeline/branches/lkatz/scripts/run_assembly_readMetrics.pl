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

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH = $FSFind::LOG || *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

exit(main());

sub main() {
  $settings = AKUtils::loadConfig($settings);
  die(usage($settings)) if @ARGV<1;

  my @cmd_options=qw(help fast qual_offset=i minLength=i);
  GetOptions($settings, @cmd_options) or die;
  $$settings{qual_offset}||=33;

  my %final_metrics;
  for my $input_file(@ARGV){
    my $file=File::Spec->rel2abs($input_file);
    die("Input or file $file not found") unless -f $file;
  
    my %metrics;
    if($$settings{fast}){
      %metrics=fastFastqStats($file,$settings);
    } else {
      %metrics=readMetrics($file,$settings);
    }
    $final_metrics{$file}=\%metrics;
  }

  printMetrics(\%final_metrics,$settings);
  return 0;
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


sub readMetrics{
  my($file,$settings)=@_;
  my $nullQual=chr($$settings{qual_offset}); # for whenever a quality score cannot be found

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
    for("$path$name$ext.qual","$path$name.qual"){
      $qualfile=$_ if(-e $_);
    }
    # convert the quality file to Illumina qual strings
    if(-e $qualfile){
      $qual=AKUtils::readMfa($qualfile,{keep_whitespace=>1});
      while(my($id,$qualstr)=each(%$qual)){
        my $newQual;
        my @tmpQual=split(/\s+/,$qualstr);
        $newQual.=chr($_+$$settings{qual_offset}) for(@tmpQual);
        $$qual{$id}=$newQual;
      }
    }
    $seqs=[values(%$seqs)];
    $qual=[values(%$qual)];
  } else {
    die "$ext extension not supported";
  }
  
  # STEP 1b: FILL IN QUALITY VALUES FOR READS WHOSE QUALITY IS NULL
  if(!@$qual){
    warn "Warning: There are no quality scores associated with $file";
    #$$qual{$_}=$nullQual x length($$seqs{$_}) for(keys(%$seqs));
    for(my $i=0;$i<@$seqs;$i++){
      $$qual[$i]=$nullQual x length($$seqs[$i]);
    }
  }

  # STEP 2: METRICS
  my $seqCounter=0;
  my $minLength=$$settings{minLength} || 0;
  my($totalReadLength,$maxReadLength,$totalReadQuality,$totalQualScores,$avgReadQualTotal)=(0,0);
  my $minReadLength=9999999999999999999999;
  #while(my($id,$seq)=each(%$seqs)){
  for(my $i=0;$i<@$seqs;$i++){
    my $seq=$$seqs[$i];
    
    # read metrics
    my $readLength=length($seq);
    next if($minLength && $readLength<$minLength);
    $totalReadLength+=$readLength;
    $maxReadLength=$readLength if($readLength>$maxReadLength);
    $minReadLength=$readLength if($readLength<$minReadLength);
    
    # quality metrics
    my $qualStr=$$qual[$i];
    my @qual;
    if(!$qualStr){
      warn "Warning: Could not find qual for seq number $i. Setting to 0.";
      $qualStr=$nullQual x $readLength;
    }
    @qual=map(ord($_)-$$settings{qual_offset},split(//,$qualStr));
    my $sumReadQuality=sum(@qual);
    my $thisReadAvgQual=$sumReadQuality/$readLength;
    $avgReadQualTotal+=$thisReadAvgQual;
    $totalReadQuality+=$sumReadQuality;
    $totalQualScores+=$readLength;
    
    $seqCounter++;
  }
  my $avgReadLength=$totalReadLength/$seqCounter;
  my $avgQuality=$totalReadQuality/$totalQualScores;
  #my $avgQualPerRead=$avgReadQualTotal/$seqCounter;

  my %metrics=(
    avgReadLength=>$avgReadLength,
    totalBases=>$totalReadLength,
    maxReadLength=>$maxReadLength,
    minReadLength=>$minReadLength,
    avgQuality=>$avgQuality,
    #avgQualPerRead=>$avgQualPerRead,
    numReads=>$seqCounter,
  );
  return %metrics;
}

sub readFastq{
  my($fastq,$settings)=@_;
  my $seqs=[];
  my $quals=[];
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

    push(@$seqs,$sequence);
    push(@$quals,$qual);
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
