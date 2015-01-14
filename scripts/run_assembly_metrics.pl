#!/usr/bin/env perl

# run_assembly_metrics.pl: finds metrics of a given assembly and prints them
# Author: Lee Katz (lskatz@gatech.edu)

package PipelineRunner;
my ($VERSION) = ('$Id: $' =~ /,v\s+(\d+\S+)/o);

my $settings = {
    appname => 'cgpipeline',
    # these are the subroutines for all assembly metrics
    metrics=>[qw(N50 genomeLength longestContig numContigs avgContigLength assemblyScore minContigLength expectedGenomeLength kmer21 GC)],
    # these are the subroutines for all standard assembly metrics
    stdMetrics=>[qw(genomeLength N50 numContigs assemblyScore)],
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
use threads;
use Thread::Queue;

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH = $FSFind::LOG || *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

exit(main());

sub main() {
  $settings = AKUtils::loadConfig($settings);
  die(usage($settings)) if @ARGV<1;

  my @cmd_options=qw(outputType=s expectedGenomeLength=i minContigLength=i statistic=s@ numberOnly allMetrics numcpus=i help);
  GetOptions($settings, @cmd_options) or die;
  $$settings{minContigLength}||=500;
  $$settings{expectedGenomeLength}||=0;
  $$settings{numcpus}||=AKUtils::getNumCPUs();
  $$settings{outputType}||="";
  die(usage($settings)) if $$settings{help};

  my @input_file;
  for my $file(@ARGV){
    $file=File::Spec->rel2abs($file);
    if(!-f $file){
      warn("Warning input or reference file $file not found") unless -f $file;
      next;
    }
    push(@input_file,$file);
  }

  if($$settings{outputType} eq 'lengths'){
    printContigDetails(\@input_file,$settings);
  } else {
    # start off the threads
    my $Q=Thread::Queue->new;
    my @thr;
    $thr[$_]=threads->new(\&assemblyMetricsWorker,$Q,$settings) for(0..$$settings{numcpus}-1);

    for my $file(@input_file){
      $Q->enqueue($file);
    }

    # close out the threads
    $Q->enqueue(undef) for(@thr);
    my %metrics;
    for(@thr){
      my $metrics=$_->join;
      %metrics=(%metrics,%$metrics);
    }

    printMetrics(\%metrics,$settings);
  }

  return 0;
}

sub assemblyMetricsWorker{
  my($Q,$settings)=@_;
  my %metrics;
  while(defined(my $file=$Q->dequeue)){
    my $metrics=assemblyMetrics($file,$settings);
    next if(keys(%$metrics) < 1);
    $metrics{$file}=$metrics;
  }
  return \%metrics;
}

sub assemblyMetrics($$){
  my($file,$settings)=@_;
  my($result);

  my $seqs={};
  if(-s $file > 0){
    $seqs=AKUtils::readMfa($file);
  }

  $seqs=filterSeqs($seqs,$settings);
  if(keys(%$seqs) < 1){
    logmsg "WARNING: Could not find any sequences in $file. Skipping";
    return {};
  }

  # put the assembly score last because it depends on everything else
  my @statsToProcess=sort({
    return 1 if($a=~/assemblyScore/i);
    return 0;
  } @{$$settings{metrics}});

  # some default metrics
  my $metrics={
    File=> $file,
    minContigLength=>$$settings{minContigLength},
    expectedGenomeLength=>$$settings{expectedGenomeLength},
  };

  # calculate metrics
  foreach my $statistic (@statsToProcess){
    my $stat;
    next if(!defined(&$statistic));
    if($statistic eq 'assemblyScore'){
      $stat=&$statistic($metrics,$settings);
    } else{
      $stat=&$statistic($seqs,$settings);
    }
    $$metrics{$statistic}=$stat;
  }
  return $metrics;

  # set the output in the result hash
  my $statsToOutput=$$settings{statistic} || $$settings{stdMetrics};
  unshift(@$statsToOutput,qw(File minContigLength expectedGenomeLength)) if(!$$settings{numberOnly});
  for my $statistic(@$statsToOutput){
    if($$settings{numberOnly}){
      $result.=$$metrics{$statistic}."\n";
    } else {
      $result.=join("\t",$statistic,$$metrics{$statistic})."\n";
    }
  }
  chomp($result);

  return $result;
}

sub printMetrics{
  my ($metrics,$settings)=@_;
  my $header=$$settings{statistic} || $$settings{stdMetrics};

  # default output
  $header=$$settings{metrics} if($$settings{allMetrics});
  unshift(@$header,"File") if(!$$settings{numberOnly});
  my $d="\t"; # field delimiter
  print join($d,@$header)."\n" if(!$$settings{numberOnly});
  while(my($file,$m)=each(%$metrics)){
    my @line;
    for(@$header){
      push(@line,$$m{$_});
    }
    print join($d,@line)."\n";
  }
  return 1;
}

sub filterSeqs{
  my($seqs,$settings)=@_;
  my $newSeqs={};
  while(my ($id,$seq)=each(%$seqs)){
    next if($$settings{minContigLength} > 0 && length($seq) < $$settings{minContigLength});
    $$newSeqs{$id}=$seq;
  }
  return $newSeqs;
}

sub printContigDetails{
  my($file,$settings)=@_;
  print join("\t",qw(file contig length GC kmer21))."\n";
  for my $f(@$file){
    my $seqs=AKUtils::readMfa($f);
    if(keys(%$seqs) < 1){
      next;
    }
    while(my($id,$seq)=each(%$seqs)){
      my $sc={$id=>$seq}; # Single Contig
      print join("\t",$f,$id,length($seq),GC($sc),kmer21($sc))."\n";
    }
  }
  return 1;
}

sub avgContigLength{
  my($seqs)=@_;
  my $length=0;
  my $numContigs=0;
  foreach my $seqname (keys %$seqs){
    $length+=length($$seqs{$seqname});
    $numContigs++;
  }
  return 0.01 if($numContigs==0);
  return int(($length/$numContigs));
}

# Find the longest contig of an assembly.
# param: seqs in Andrey's format
# return int The contig size, in bp
sub longestContig($){
  my($seqs)=@_;
  my $longest=0;
  foreach my $seqname (keys %$seqs){
    my $contigLength=length($$seqs{$seqname});
    if($contigLength>$longest){
      $longest=$contigLength;
    }
  }
  $longest||=0.01;
  return $longest;
}

# Find the N50 of an assembly.  The N50 is the size N such that 50% of the genome is contained in contigs of size N or greater.
# param: fasta filename
# optional param: $$settings
# return int The N50 in bp.
sub N50($;$){
  my($seqs,$settings)=@_;
  my($seqname,@seqs,$numSeqs,$N50);
  # put the seqs into an array. No defline needed in this subroutine, so it can be discarded
  foreach $seqname (keys %$seqs){
    push(@seqs,$$seqs{$seqname});
  }
  $numSeqs=@seqs;
  # order the contigs by size
  # Largest contigs are towards $seqs[0], and smallest contigs are at the last elements
  @seqs=sort{
    length($a)<=>length($b);
  } @seqs;
  # if the size is not provided, assume the size of the assembly is the genome size
  my $genomeSize = $$settings{expectedGenomeLength}||genomeLength($seqs);
  my $halfGenomeSize=$genomeSize/2;
  my $currentGenomeLength=0;
  # find the contig that resides at (size/2)
  for(my $i=0;$i<$numSeqs;$i++){
    my $seqLength=length($seqs[$i]);
    $currentGenomeLength+=$seqLength;
    if($currentGenomeLength>$halfGenomeSize){
      $N50=$seqLength;
      last;
    }
  }
  # return the length of the contig
  $N50||=0.01;
  return $N50;
}

# what percent of kmers are repeated
sub kmer21{
  my($seqs,$settings)=@_;
  my %mer21;
  for my $sequence(values(%$seqs)){
    my $length=length($sequence);
    for(my $i=0;$i<$length-21;$i++){
      my $kmer=substr($sequence,$i,21);
      $mer21{$kmer}++;
    }
  }

  my $mer21;
  my $totalKmer;
  for(values(%mer21)){
    $mer21++ if($_>1);
    $totalKmer+=$_;
  }
  return 0 if ($mer21<=0);
  my $kmer21Freq=sprintf("%0.4f",($mer21/$totalKmer));
  return $kmer21Freq;
}

# Find the length of a genome assembly
# param: fasta filename
# return int The size of the genome in bp
sub genomeLength($){
  my($seqs)=@_;
  my ($seqname);
  my $length=0;
  foreach $seqname (keys %$seqs){
    $length+=length($$seqs{$seqname});
  }
  $length||=0.01;
  return $length;
}

# percent GC
sub GC($){
  my($seqs)=@_;
  my($gc,$length);
  foreach my $seqname(keys %$seqs){
    $length+=length($$seqs{$seqname});

    my $sequence=$$seqs{$seqname};
    my $s=$sequence; # I think this should help avoid modifying the sequence
    $gc+=($s=~s/[GCgc]//g);
  }
  my $GC=$gc/$length;
  $GC=(sprintf("%0.3f",$GC) * 100) . '%';
  return $GC;
}

# make a crude assembly score for comparing the best assemblies
sub assemblyScore{
  my($metrics,$settings)=@_;
  my $percentGenomeCovered=1;
  if($$settings{expectedGenomeLength}){
    $percentGenomeCovered=1-abs(1-$$metrics{genomeLength}/$$settings{expectedGenomeLength});
  }
  $percentGenomeCovered=1/9999999999 if($percentGenomeCovered<=0);
  my $score=$$metrics{N50}/$$metrics{numContigs} * $percentGenomeCovered;
  die "score is negative and so it cannot be put on the log scale because: $$metrics{N50}/$$metrics{numContigs} * [1-abs(1-$$settings{expectedGenomeLength}/$$metrics{genomeLength})]\n" if($score<=0);
  my $logscore=log($score);
  $logscore=sprintf("%0.3f",$logscore); # round to 1000s place
  return $logscore;
}


# Find the number of contigs in an assembly
# param: fasta filename
# return int The number of contigs
sub numContigs($){
  my($seqs)=@_;
  my $numContigs=scalar(values(%$seqs));
  $numContigs=9999999999 if($numContigs<1);
  
  return $numContigs;
}

sub usage{
  my ($settings)=@_;
  my $text="Prints useful assembly statistics
  Usage: $0 assembly1.fasta [-e expectedGenomeLength] [-m minContigLength]
  -e genome length Helps with N50 calculation and assemblyScore but not necessary
  -m size   Only consider contigs of size m in the calculations
  -number   Output only the number of the metric(s) you supplied and not the header. Works well with -s
  --numcpus 1 The number of threads to use
  -s stat   (in contrast, see -a)
    only print out the value for this stat (or these stats) only.
    Possible values: ".join(", ",@{$$settings{metrics}})."
    Default output: ".join(", ",@{$$settings{stdMetrics}})."
  -a to output all stats. To select only some stats, use -s
  -h for more help";
  return $text if(!$$settings{help});
  $text.="\n  ADDITIONAL HELP
    -o outputType Make different output type.  Output type can be one of the following
      details - metrics on individual contigs
      default - no change
    NOTES
    The assembly score is calculated as a log of (N50/numContigs * percentOfGenomeCovered)
    The percent of the genome covered is 1 if you do not supply an expected genome length.
    The percent of the genome covered counts against you if you assemble higher than the expected genome size.
  EXAMPLES - combine piped commands like column sort to get even better output
    $0 *.fasta | column -t             # pretty output
    $0 *.fasta | sort -k5,5n           # sorted output
    $0 *.fasta -a                      # get all stats
    $0 *.fasta | sed 's|.*/||'         # remove directories from filenames
    N50=`$0 file.fasta -number -s N50` # set a bash variable equal to an assembly's N50
    # histogram of contig lengths
    $0 file.fasta -o details | perl -lane 'print(int(\$F[2]/100000)*100000)'|sort -n|uniq -c
  ";
  return $text;
}

