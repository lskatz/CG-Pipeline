#!/usr/bin/env perl

# run_assembly_trimClean: trim and clean a set of raw reads
# Author: Lee Katz <lkatz@cdc.gov>

package PipelineRunner;
my ($VERSION) = ('$Id: $' =~ /,v\s+(\d+\S+)/o);

my $stats;

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
use POSIX;

use Data::Dumper;
use threads;
use Thread::Queue;
use threads::shared;

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {print STDERR "$0: ".(caller(1))[3].": @_\n";}

my %threadStatus:shared;
my @IOextensions=qw(.fastq .fastq.gz .fq);

exit(main());

sub main{
  logmsg "WARNING: Removing duplicates is experimental";
  my $settings = {
      appname => 'cgpipeline',
  };

  # get settings from cg_pipeline/conf/cgpipelinerc and ./cgpipelinerc
  $settings=AKUtils::loadConfig($settings);
  # get CLI flags into $settings with Getopt::Long
  GetOptions($settings,qw(help tempdir=s)); 
  die usage() if(@ARGV<1 || $$settings{help});
  
  # additional settings using ||= operator
  $$settings{tempdir}||=AKUtils::mktempdir(); # any temp file you make should go under the tempdir
  $$settings{is_PE}||=AKUtils::is_fastqPE($infile,$settings); # 0 for SE; 1 for PE
  $$settings{poly}=$$settings{is_PE}+1;
  
  my($peFastq,$seFastq)=removeDuplicateReads($infile,$settings);

  return 0;
}

sub removeDuplicateReads{
  my($infile,$settings)=@_;

  my @nodupes;
  for my $fastq(@$infile){
    my ($basename,$dir,$inSuffix) = fileparse($fastq,@IOextensions);
    my $poly=$$settings{poly}||checkPolyness($fastq,$settings);
    my $nodupes="$$settings{tempdir}/$basename.noduplicates.fastq";
    push(@nodupes,$nodupes);
    logmsg "$fastq => $nodupes";
    open(OUT,">",$nodupes) or die "Could not open fastq file $nodupes for writing: $!";
    if($inSuffix=~/\.fastq$/){
      open(FQ,'<',$fastq) or die "Could not open $fastq for reading: $!";
    }
    elsif($inSuffix=~/\.fastq\.gz/){
      open(FQ,"gunzip -c $fastq | ") or die "Could not open $fastq for reading: $!";
    }
    else{
      die "Could not determine the file type for reading based on your extension $inSuffix";
    }
    my $i=0;
    my %read;
    while(my $id=<FQ>){
      my $sequence=<FQ>; chomp($sequence);
      my $discard=<FQ>;
      my $qual=<FQ>;
      my $read="$id$sequence\n+\n$qual";
      my $hashId=$sequence;
      for(my $j=1;$j<$poly;$j++){
        my $id2=<FQ>;
        my $sequence2=<FQ>; chomp($sequence2);
        my $discard2=<FQ>;
        my $qual2=<FQ>;
        $read.="$id2$sequence2\n+\n$qual2";
        $hashId.="~~~~$sequence2"; # tildes are probably not in the sequence, so it's probably safe for putting into the hash ID
      }
      $read{$hashId}=$read;
    }
    close FQ;

    while(my($hashId,$read)=each(%read)){
      print OUT $read;
    }
    close OUT;
  }
  logmsg "Done";
  return \@nodupes;
}


sub usage{
  "Removes duplicate reads from a raw read set
  Usage: $0 read.fastq[.gz] > read.fastq.gz
  ";
}

