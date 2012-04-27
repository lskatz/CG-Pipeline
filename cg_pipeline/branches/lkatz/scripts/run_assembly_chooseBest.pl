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
  die(usage($settings)) if @ARGV<1;

  my @cmd_options=qw(output=s expectedGenomeLength=i metrics);
  GetOptions($settings, @cmd_options) or die;

  $$settings{outfile} = $$settings{output} || "$0.out.fasta";
  $$settings{outfile} = File::Spec->rel2abs($$settings{outfile});
  $$settings{min_contig_length}||=1;
  
  my @input_files=@ARGV;
  logmsg "Determining best assembly from ".join(", ",@input_files).".";
  foreach my $file (@input_files) {
	  $file = File::Spec->rel2abs($file);
	  die("Input or reference file $file not found") unless -f $file;
  }
  my $final_seqs=bestAssemblySeqs(\@input_files,$settings);
  
  AKUtils::printSeqsToFile($final_seqs, $$settings{outfile}, {order_seqs_by_name => 1});
  logmsg "Output is in $$settings{outfile}";

  if($$settings{metrics}){
    logmsg "Gathering metrics...";
    my $metricsCommand="run_assembly_metrics.pl $$settings{outfile} -m $$settings{min_contig_length}";
    $metricsCommand.=" -e $$settings{expectedGenomeLength}" if($$settings{expectedGenomeLength});
    system($metricsCommand);
  }

  return 0;
}

# Find the best assembly, given the contigs in that assembly
# params: list each assembly, composed of contigs
# return seqs of best assembly
sub bestAssemblySeqs($$){
  my($seqs,$settings)=@_;
  my($i,%metrics);

  for(@$seqs){
    $metrics{$_}=`run_assembly_metrics.pl '$_' -n -s assemblyScore`+0;
  }
  my @metrics=values(%metrics);
  my @files=keys(%metrics);
  my $bestScore=shift(@metrics);
  my $bestFile=shift(@files);
  for($i=0;$i<@files;$i++){
    next if($bestScore>$metrics[$i]);
    $bestScore=$metrics[$i];
    $bestFile=$files[$i];
  }
  logmsg "Best assembly is $bestFile with a score of $bestScore";
  return AKUtils::readMfa($bestFile);
}

sub usage{
  "Usage: $0 assembly1.fasta [,assembly2.fasta...] -output bestAssembly.fasta [-e expectedGenomeLength]
  -e expectedGenomeLength in bp
  -m
    if you want final metrics to be displayed
  Tip: put your preffered assembly first to help break ties.
  "
}
