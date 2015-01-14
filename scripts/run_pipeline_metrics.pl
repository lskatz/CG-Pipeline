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
  die("  Usage: $0 projectName [projectName2...]\n") if (@ARGV<1);

  my @cmd_options=('output=s','project=s@');
  GetOptions($settings, @cmd_options) or die;

  my $project=$$settings{project} || ();
  push(@$project,@ARGV); # pick up any remaining projects that were given as ARGV
  for my $p(@$project){
    $p=File::Spec->rel2abs($p);
    die("Project $project not found") unless -d $p;
    metrics($p,$settings);
  }
  return 0;
}

sub metrics{
  my($project,$settings)=@_;
  
  if(-e "$project/assembly.fasta"){
    system("run_assembly_metrics.pl $project/assembly.fasta");
    die if $?;
  }
  if(-e "$project/prediction.gb"){
    system("run_prediction_metrics.pl $project/prediction.gb");
    die if $?;
  }
  if(-e "$project/annotation.gb"){
    system("run_prediction_metrics.pl $project/annotation.gb");
    die if $?;
  }
}

