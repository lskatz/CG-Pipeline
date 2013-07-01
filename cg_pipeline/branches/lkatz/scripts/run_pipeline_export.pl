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
use Archive::Tar;
use File::Slurp;

sub logmsg {my $FH = *STDERR; print $FH "$0: ".(caller(1))[3].": @_\n";}
exit(main());

sub main{
  die usage() if @ARGV<1;
  GetOptions($settings,qw(project=s@ outfile=s));

  my $project=$$settings{project} or die "Error: no project given\n".usage();
  
  $$settings{tempdir}||=AKUtils::mktempdir();
  my @filesToCompress;
  for my $p (@$project){
    my $file=locateFile($p,$settings);
    my($name,$path,$suffix)=fileparse($file,qw(.gb .gbk .fasta));
    $p=~s/\/+$//; # remove trailing slashes
    my $tmpfile="$$settings{tempdir}/$p$suffix";
    copy($file,$tmpfile) or die "Could not copy $file to $tmpfile because $!\n";
    push(@filesToCompress,$tmpfile);
  }

  if(my $outfile=$$settings{outfile}){
    my $tar=Archive::Tar->new;
    for my $file(@filesToCompress){
      my($name,$path,$suffix)=fileparse($file,qw(.gb .gbk .fasta));
      my $fileContents=read_file($file);
      logmsg "Adding cgpExport/$name$suffix to the archive";
      $tar->add_data("cgpExport/$name$suffix",$fileContents);
      die "Could not add data to tar archive" if $?;
    }
    $tar->write($outfile,COMPRESS_GZIP);
    die "Error: Problem with 'tar zcvf'" if $?;
  } else {
    die "STDOUT not supported right now";
  }
  #printFile($file,$settings);

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
  die "Could not locate a genome project file in $project out of ".join(" ",@desiredFile);
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
  Usage: $0 -p project [-p project2 ...] -o outfile.tar.gz
  -p a cg_pipeline project directory
  -o outfile
    The file will be compressed into a tar.gz file
    Uncompress with 'tar zxvf file.tar.gz' and the results will be in the cgpExport directory
  ";
}
