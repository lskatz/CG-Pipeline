#!/usr/bin/env perl

# run_pipeline_metrics.pl: finds metrics for a pipeline project. 
# Is a wrapper around run_*_metrics.pl 
# Author: Lee Katz ()

package PipelineRunner;
my ($VERSION) = ('$Id: $' =~ /,v\s+(\d+\S+)/o);

my $settings = {
    appname => 'cgpipeline',
};
my $stats;

use strict;
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
  die("  Usage: $0 -p projectName\n") if @ARGV<1;

  my @cmd_options=('output=s','project=s');
  GetOptions($settings, @cmd_options) or die;

  my $project=$$settings{project};
  my $project=File::Spec->rel2abs($project);
  die("Project $project not found") unless -d $project;

  my $metrics=metrics($project,$settings);
  # TODO print metrics in an understandable and parsable way
  print join("\t","Project",$project)."\n";
  foreach my $metric (sort keys %$metrics){
    print join("\t",$metric,$$metrics{$metric});
    print "\n";
  }
  return 0;
}

sub metrics{
  my($project,$settings)=@_;
  my (%metrics,$command,$out);
  
  return \%metrics if(!-e "$$settings{project}/assembly.fasta");
  $out=`run_assembly_metrics.pl $$settings{project}/assembly.fasta`;
  readMetrics($out,\%metrics,'assembly_');
  
  return \%metrics if(!-e "$$settings{project}/prediction.gb");
  $out=`run_prediction_metrics.pl $$settings{project}/prediction.gb`;
  readMetrics($out,\%metrics,'prediction_');

  return \%metrics if(!-e "$$settings{project}/annotation.gb");
  $out=`run_prediction_metrics.pl $$settings{project}/annotation.gb`;
  readMetrics($out,\%metrics,'annotation_');

  return \%metrics;
}

# read a metrics output into a hash
sub readMetrics{
  my($metricsStr,$hash,$key_prefix)=@_;
  
  for(split /\n/,$metricsStr){
    my($key,$value)=split /\t/;
    $$hash{"$key_prefix$key"}=$value;
  }
  return 1;
}
