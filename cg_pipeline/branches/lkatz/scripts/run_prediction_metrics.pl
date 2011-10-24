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
sub logmsg {my $FH = $FSFind::LOG || *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

exit(main());

sub main() {
  $settings = AKUtils::loadConfig($settings);
  die("  Usage: $0 prediction.gb\n") if @ARGV<1;

  my @cmd_options=('output=s');
  GetOptions($settings, @cmd_options) or die;

  my $input_file=$ARGV[0];
  my $file=File::Spec->rel2abs($input_file);
  die("Input or reference file $file not found") unless -f $file;

  my $metrics=predictionMetrics($file,$settings);
  # TODO print metrics in an understandable and parsable way
  print join("\t","File",$input_file)."\n";
  foreach my $metric (keys %$metrics){
    print join("\t",$metric,$$metrics{$metric})."\n";
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
      #for my $tag ($feat->get_all_tags){
      #  for my $value($feat->get_tag_values($tag)){
      #  }
      #}
    }
  }
  return $metrics;
}

