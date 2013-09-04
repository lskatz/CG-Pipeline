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

$0=fileparse $0;
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; print STDERR ("$0: ".(caller(1))[3].": ".$e); exit -1; };
sub logmsg {return 0 if $$settings{quiet}; my $FH = *STDERR; print $FH "$0: ".(caller(1))[3].": @_\n";}
exit(main());

sub main{
  die usage() if @ARGV<1;
  GetOptions($settings,qw(export clean casm cpred cann boolean quiet));
  die usage() if @ARGV<1;

  my @project=@ARGV;

  logmsg "WARNING: exit codes might not work as intended when you get to 127 genome projects and above" if(@project>=127);

  my $exit_code=0;
  for my $p(@project){
    if(!is_cgp_project($p,$settings)){
      logmsg "$p is not a CGP project. I will skip it.";
      next;
    } else { 
      logmsg "CGP project $p";
      $exit_code++;
    }
    exportProject($p,$settings) if($$settings{export});
    cleanProject($p,$settings) if($$settings{clean});
    cleanAssembly($p,$settings) if($$settings{casm});
    cleanPrediction($p,$settings) if($$settings{cpred});
    cleanAnnotation($p,$settings) if($$settings{cann});
  }
  return $exit_code;
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
  my $return=0;
  $return+=cleanAssembly($project,$settings); 
  $return+=cleanPrediction($project,$settings); 
  $return+=cleanAnnotation($project,$settings);
  return $return;
}

sub cleanAssembly{
  my($project,$settings)=@_;
  my $dir="$project/build/assembly";
  return 0 if(!glob("$dir/*"));
  File::Path->remove_tree($dir,{keep_root=>1,error=>\my $err});
  if(@$err){
    logmsg "At least one error was returned when cleaning. Attempting one more time.";
    File::Path->remove_tree($dir,{keep_root=>1,error=>$err});
    if(@$err){
      logmsg "I still encountered at least one error.";
      die Dumper $err;
    }
  }
  return 1;
}
sub cleanPrediction{
  my($project,$settings)=@_;
  my $dir="$project/build/prediction";
  return 0 if(!glob("$dir/*"));
  File::Path->remove_tree($dir,{keep_root=>1,error=>\my $err});
  if(@$err){
    logmsg "At least one error was returned when cleaning. Attempting one more time.";
    File::Path->remove_tree($dir,{keep_root=>1,error=>$err});
    if(@$err){
      logmsg "I still encountered at least one error.";
      die Dumper $err;
    }
  }
  return 1;
}
sub cleanAnnotation{
  my($project,$settings)=@_;
  my $dir="$project/build/annotation";
  return 0 if(!glob("$dir/*"));
  File::Path->remove_tree($dir,{keep_root=>1,error=>\my $err});
  if(@$err){
    logmsg "At least one error was returned when cleaning. Attempting one more time.";
    File::Path->remove_tree($dir,{keep_root=>1,error=>$err});
    if(@$err){
      logmsg "I still encountered at least one error.";
      die Dumper $err;
    }
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
  "Manages a CG-Pipeline project. With no options, simply reports whether a directory is a CGP project.
  Usage: $0 project [project2 ...] -e > outfile
  project: a cg_pipeline project directory (created by 'run_pipeline create')
  -e export the project directory
  --clean clean out all temporary files, 
          or use one or more of the following to clean a specific stage
    -casm -cpre -cann
  -q for quiet

  Exit codes can be used to your advantage with this script
  for example:
    system('$0 -q CGPdir'); \$e=\$? >> 8; if(\$e){ print 'This is a CGP project directory!';}
  However, an exit code of 255 signifies an error.
  ";
}
