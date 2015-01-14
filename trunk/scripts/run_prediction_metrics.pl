#!/usr/bin/env perl

# run_prediction_metrics.pl: finds metrics of a given prediction genbank file
# Author: Lee Katz (lskatz@gatech.edu)

package PipelineRunner;
my ($VERSION) = ('$Id: $' =~ /,v\s+(\d+\S+)/o);

my $settings = {
    appname => 'cgpipeline',
};
my $stats;

use strict;
no strict "refs";
use warnings;
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
use Bio::Perl;

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };

exit(&main());

sub main() {
  $settings = AKUtils::loadConfig($settings);
  die("  Usage: $0 prediction.gb [prediction2.gb...]\n") if @ARGV<1;

  my @cmd_options=('output=s');
  GetOptions($settings, @cmd_options) or die;

  my %header;
  my %output;
  my @input_file=@ARGV;
  for my $input_file(@input_file){
    my $file=File::Spec->rel2abs($input_file);
    die("Input or reference file $file not found") unless -f $file;

    my $metrics=predictionMetrics($file,$settings);
    $output{$file}=$metrics;
    @header{keys(%$metrics)}=keys(%$metrics); # to get a unique set of headers
    $input_file=$file; # change this by reference, for later
  }

  # print results using headers found while getting the metrics
  delete($header{File});
  my @header=keys(%header);
  unshift(@header,"File"); # make sure File is first

  # print the header first, then the metrics for each prediction file
  print join("\t",@header)."\n";
  for my $input_file(@input_file){
    for my $h(@header){
      my $value=$output{$input_file}{$h} || ".";
      print "$value\t";
    }
    print "\n";
  }
  return 0;
}

sub predictionMetrics($$){
  my($file,$settings)=@_;
  my $metrics; # hash of prediction metrics

  my %seenSeq;
  my $seqio=Bio::SeqIO->new(-file=>$file,-format=>'genbank');
  while(my $seq=$seqio->next_seq){
    next if($seenSeq{$seq->id}++);
    $$metrics{genomeLength}+=$seq->length;
    for my $feat ($seq->get_SeqFeatures){
      $$metrics{$feat->primary_tag}++;
    }
  }
  $$metrics{File}=$file;
  return $metrics;
}

