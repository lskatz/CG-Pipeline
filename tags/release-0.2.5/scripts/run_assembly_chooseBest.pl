#!/usr/bin/env perl

# run_assembly_chooseBest: Choose the best assembly based on several metrics derived from their FASTA files
# Author: Lee Katz (lskatz@gatech.edu)

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
  die("  Usage: $0 assembly1.fasta [,assembly2.fasta...] -output bestAssembly.fasta\n  Tip: put your preffered assembly first to help break ties.\n") if @ARGV<1;

  my @cmd_options=('output=s');
  GetOptions($settings, @cmd_options) or die;

  $$settings{outfile} = $$settings{output} || "$0.out.fasta";
  $$settings{outfile} = File::Spec->rel2abs($$settings{outfile});
  
  my @input_files=@ARGV;
  logmsg "Determining best assembly from ".join(", ",@input_files).".";
  foreach my $file (@input_files) {
	  $file = File::Spec->rel2abs($file);
	  die("Input or reference file $file not found") unless -f $file;
  }
  my $final_seqs=bestAssemblySeqs(\@input_files);
  
  AKUtils::printSeqsToFile($final_seqs, $$settings{outfile}, {order_seqs_by_name => 1});
  logmsg "Output is in $$settings{outfile}.";

  return 0;
}

# Find the best assembly, given the contigs in that assembly
# params: list each assembly, composed of contigs
# return seqs of best assembly
sub bestAssemblySeqs($){
  my($seqs)=@_;
  my(@stats,@vote,$winningSeqIndex,$i,@seqsFilename);

  for($i=0;$i<@$seqs;$i++){
	$seqsFilename[$i]=$$seqs[$i];
    $$seqs[$i]=AKUtils::readMfa($$seqs[$i]);
  }

  # for each statistic give a point to the winner (BIGGEST number)
  foreach my $statistic qw(N50 genomeLength longestContig){
    my $largestStat=-1000; # this can be any very negative number that all assemblies will beat
    $winningSeqIndex=0; # by default, the first assembly wins the vote
    for($i=0;$i<@$seqs;$i++){
      my $stat=&$statistic($$seqs[$i]);
      if($stat>$largestStat){
        $winningSeqIndex=$i;
        $largestStat=$stat;
      }
    }
    logmsg "$seqsFilename[$winningSeqIndex] wins a point with metric $statistic ($largestStat).";
    $vote[$winningSeqIndex]++;
  }
  # for each statistic give a point to the winner (LOWEST number)
  foreach my $statistic qw(numContigs){
    my $smallestStat=100000000000000; # this can be any very positive number that all assemblies will beat
    $winningSeqIndex=0; # by default, the first assembly wins the vote
    for($i=0;$i<@$seqs;$i++){
      my $stat=&$statistic($$seqs[$i]);
      if($stat<$smallestStat){
        $winningSeqIndex=$i;
        $smallestStat=$stat;
      }
    }
    logmsg "$seqsFilename[$winningSeqIndex] wins a point with metric $statistic ($smallestStat).";
    $vote[$winningSeqIndex]++;
  }

  # Which has the most votes?  The most votes is the best assembly
  # TODO  In the case of a tie for first place, run this subroutine one more time with the top two winners.
  # TODO  In the case where the tiebreaker didn't work, introduce another metric
  my $mostVotes=-1000;
  $winningSeqIndex=0; # by default, the first assembly wins the vote
  for(my $i=0;$i<@vote;$i++){
    if($vote[$i]>$mostVotes){
      $mostVotes=$vote[$i];
      $winningSeqIndex=$i;
    }
  }
  logmsg "The best assembly is $seqsFilename[$winningSeqIndex] with $vote[$winningSeqIndex] votes.";

  return $$seqs[$winningSeqIndex];
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
  return $longest;
}

# Find the N50 of an assembly.  The N50 is the size N such that 50% of the genome is contained in contigs of size N or greater.
# param: fasta filename
# optional param: total size of genome
# return int The N50 in bp.
sub N50($;$){
  my($seqs,$genomeSize)=@_;
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
  $genomeSize ||= genomeLength($seqs);
  my $halfGenomeSize=$genomeSize/2;
  my $currentGenomeLength=0;
  # find the contig that resides at (size/2)
  for(my $i=0;$i<$numSeqs;$i++){
    my $seqLength=length($seqs[$i]);
    $currentGenomeLength+=$seqLength;
#   print join ("\t", ($i,$seqLength,$currentGenomeLength,$halfGenomeSize))."\n";
    if($currentGenomeLength>$halfGenomeSize){
      $N50=$seqLength;
      last;
    }
  }
  # return the length of the contig
  return $N50;
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
  return $length;
}

# Find the number of contigs in an assembly
# param: fasta filename
# return int The number of contigs
sub numContigs($){
  my($seqs)=@_;
  my(@seqValue)=values %$seqs;
  return scalar @seqValue;
}

