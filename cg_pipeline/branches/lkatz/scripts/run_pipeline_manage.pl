#!/usr/bin/env perl

# run_pipeline_manage: manage a CGP project directory
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
use Archive::Tar;
use File::Slurp;

sub logmsg {my $FH = *STDERR; print $FH "$0: ".(caller(1))[3].": @_\n";}
exit(main());

sub main{
  die usage() if @ARGV<1;
  GetOptions($settings,qw(export clean));
  die usage() if @ARGV<1;

  my @project=@ARGV;

  for my $p(@project){
    die "ERROR: $p is not a CGP project" if(!is_cgp_project($p,$settings));
    exportProject($p,$settings) if($$settings{export});
    cleanProject($p,$settings) if($$settings{clean});
  }
  return 0;
}

sub is_cgp_project{
  my($p,$settings)=@_;
  for("$p","$p/annotation","$p/log","$p/build","$p/build/assembly","$p/build/prediction","$p/build/annotation"){
    return 0 if(!-d $_);
  }
  return 1;
}

sub cleanProject{
  my($project,$settings)=@_;
  for(glob("$project/build/assembly/* $project/build/prediction/* $project/build/annotation/*")){
    File::Path->remove_tree($_,{keep_root=>1}) or die "Could not remove $_: $!";
  }
  return 1;
}
  
sub exportProject{
  my($project,$settings)=@_;
  my $file=locateFile($project,$settings);
  my $fileContents=read_file($file);
  print $fileContents;
  return length($fileContents);
}

# locates the most advanced file: assembly, prediction, or annotation
sub locateFile{
  my($project,$settings)=@_;
  my @desiredFile=qw(annotation.gb prediction.gb assembly.fasta);
  for(@desiredFile){
    my $file="$project/$_";
    if(-e $file && -s $file>0){
      return $file;
    }
  }
  die "Could not locate a genome project file in $project out of ".join(" ",@desiredFile);
}

sub usage{
  "Manages a genome annotation project (created by 'run_pipeline create')
  Usage: $0 project [project2 ...] -e > outfile
  project: a cg_pipeline project directory
  -e export the project directory
  -c clean the project directory
  ";
}
