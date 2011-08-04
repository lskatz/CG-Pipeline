#!/usr/bin/env perl

# run_assembly_metrics.pl: finds metrics of a given assembly and prints them
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
  die("  Usage: $0 assembly1.fasta [-e expectedGenomeLength]\n") if @ARGV<1;

  my @cmd_options=('output=s','expectedGenomeLength=i');
  GetOptions($settings, @cmd_options) or die;

  my $input_file=$ARGV[0];
  my $file=File::Spec->rel2abs($input_file);
  die("Input or reference file $file not found") unless -f $file;

  my $metrics=assemblyMetrics($file,$settings);
  # TODO print metrics in an understandable and parsable way
  print "File\t$file\n";
  print "ExpectedGenomeLength\t$$settings{expectedGenomeLength}\n" if($$settings{expectedGenomeLength});
  print "$metrics\n";

  return 0;
}

sub assemblyMetrics($$){
  my($seqs,$settings)=@_;
  my($result);
  $seqs=AKUtils::readMfa($seqs);

  foreach my $statistic qw(N50 genomeLength longestContig numContigs){
    my $stat=&$statistic($seqs,$settings);
    $result.=join("\t",$statistic,$stat)."\n";
  }
  chomp($result);
  return $result;
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

