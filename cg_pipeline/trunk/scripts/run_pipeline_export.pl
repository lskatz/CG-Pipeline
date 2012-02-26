#!/usr/bin/env perl

# run_pipeline_export: export results from the pipeline
# Author: Lee Katz <lkatz@cdc.gov>

package PipelineRunner;
my ($VERSION) = ('$Id: $' =~ /,v\s+(\d+\S+)/o);

my $settings = {
  appname => 'cgpipeline',
};

use strict;
use FindBin;
use lib "$FindBin::RealBin/../lib";
$ENV{PATH} = "$FindBin::RealBin:".$ENV{PATH};
use AKUtils;

use Getopt::Long;
use File::Temp ('tempdir');
use File::Path;
use File::Spec;
use File::Copy;
use File::Basename;
use List::Util qw(min max sum shuffle);
use CGPipelineUtils;
use Data::Dumper;

sub logmsg {my $FH = *STDERR; print $FH "$0: ".(caller(1))[3].": @_\n";}
exit(main());

sub main{
  die usage() if @ARGV<1;
  GetOptions($settings,qw(project=s));

  my $project=$$settings{project} or die "Error: no project given\n".usage();
  my $file=locateFile($project,$settings);
  printFile($file,$settings);

  return 0;
}

sub locateFile{
  my($project,$settings)=@_;
  my @desiredFile=qw(annotation.gb prediction.gb assembly.fasta);
  for(@desiredFile){
    my $file="$project/$_";
    if(-e $file && -s $file>0){
      return $file;
    }
  }
  die "Could not locate a genome project file out of ".join(" ",@desiredFile);
}

sub printFile{
  my($file,$settings)=@_;
  open(IN,$file) or die "Error: Could not open file $file: $!";
  while(my $line=<IN>){
    print $line;
  }
  close IN;
  return 1;
}

sub usage{
  "Exports the results from a genome annotation project (created by 'run_pipeline create')
  Usage: $0 -p project > file
  ";
}
