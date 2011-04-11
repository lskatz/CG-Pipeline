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

  my @cmd_options=('output=s','expectedGenomeLength=i');
  GetOptions($settings, @cmd_options) or die;

  $$settings{outfile} = $$settings{output} || "$0.out.fasta";
  $$settings{outfile} = File::Spec->rel2abs($$settings{outfile});
  
  my @input_files=@ARGV;
  my @tmpFiles=();
  for(my $i=0;$i<@input_files;$i++){
    if(-s $input_files[$i]<1){
      logmsg "I will not consider the empty file $input_files[$i]";
    }
    else{
      push(@tmpFiles,$input_files[$i]);
    }
  }
  @input_files=@tmpFiles;
  logmsg "Determining best assembly from ".join(", ",@input_files).".";
  foreach my $file (@input_files) {
	  $file = File::Spec->rel2abs($file);
	  die("Input or reference file $file not found") unless -f $file;
  }
  my $final_seqs=bestAssemblySeqs(\@input_files,$settings);
  
  AKUtils::printSeqsToFile($final_seqs, $$settings{outfile}, {order_seqs_by_name => 1});
  logmsg "Output is in $$settings{outfile}";

  return 0;
}

# Find the best assembly, given the contigs in that assembly
# params: list each assembly, composed of contigs
# return seqs of best assembly
# TODO if one of the metrics is an order of magnitude higher than any other assembly's metric then give it a +1
#   This would give a distict advantage to much better assemblies (watch out for numContigs=0 and the like)
sub bestAssemblySeqs($$){
  my($seqs,$settings)=@_;
  my($i,%metrics,%vote);
  my @largestStats=qw(N50 genomeLength longestContig);
  my @smallestStats=qw(numContigs);

  # generate the metrics
  @$seqs=sort{$a cmp $b} @$seqs;
  for($i=0;$i<@$seqs;$i++){
    my %seqMetric;
    my $command="run_assembly_metrics.pl $$seqs[$i]";
    $command.=" -e $$settings{expectedGenomeSize}" if($$settings{expectedGenomeSize});
    logmsg "COMMAND $command";
    my $seqMetric=`$command`;
    my @seqMetric=split(/\n/,$seqMetric);
    for my $line (@seqMetric){
      my($key,$value)=split(/\t/,$line);
      $seqMetric{$key}=$value;
    }
    $metrics{$seqMetric{File}}=\%seqMetric;
  }

  # largest values give a vote
  foreach my $statistic (@largestStats,@smallestStats){
    my $key="---";
    # make a hash of this statistic
    my %hash;
    while(my ($assemblyFilename,$stats)=each(%metrics)){
      $hash{$assemblyFilename}=$$stats{$statistic};
    }
    $key=maxValueIndex(\%hash) if(in_array($statistic,\@largestStats));
    $key=minValueIndex(\%hash) if(in_array($statistic,\@smallestStats));
    $vote{$key}++;
  }

  my $bestAssemblyFilename=winner(\%vote);

  return AKUtils::readMfa($bestAssemblyFilename);
}

# returns the key with the biggest value
# There's probably a more optimal way to do this but it's not necessary especially since the number of assemblies will in the neighborhood of 1 to 5.  Even 100 assemblies would probably be negligable.
sub winner{
  my ($hash)=@_;
  my $mostVotes=0;
  my $winningKey="----";
  while( my($key,$value)=each(%$hash)){
    if($value>$mostVotes){
      $mostVotes=$value;
      $winningKey=$key;
    }
  }
  return $winningKey;
}
# returns the index of the array with the lowest numerical value
sub minValueIndex{
  my($hash)=@_;
  my $lowestValue=999999999; # some really large number
  my $lowestKey=-1;
  while( my($key,$value)=each(%$hash)){
    if($value<$lowestValue){
      $lowestValue=$value;
      $lowestKey=$key;
    }
  }
  return $lowestKey;
}
# largest value of an array
sub maxValueIndex{
  my($hash)=@_;
  my $largestValue=-999999999; # some really negative number
  my $largestKey=-1;
  while( my($key,$value)=each(%$hash)){
    if($value>$largestValue){
      $largestValue=$value;
      $largestKey=$key;
    }
  }
  return $largestKey;
}

# http://www.go4expert.com/forums/showthread.php?t=8978
sub in_array {
     my ($search_for,$arr) = @_;
     my %items = map {$_ => 1} @$arr; # create a hash out of the array values
     return (exists($items{$search_for}))?1:0;
}

