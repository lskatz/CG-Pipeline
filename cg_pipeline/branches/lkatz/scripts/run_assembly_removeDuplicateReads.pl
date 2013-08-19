#!/usr/bin/env perl 

# run_assembly_trimClean: trim and clean a set of raw reads
# Author: Lee Katz <lkatz@cdc.gov>

package PipelineRunner;
my ($VERSION) = ('$Id: $' =~ /,v\s+(\d+\S+)/o);

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
  my $settings = {
      appname => 'cgpipeline',
  };

  # get settings from cg_pipeline/conf/cgpipelinerc and ./cgpipelinerc
  $settings=AKUtils::loadConfig($settings);
  # get CLI flags into $settings with Getopt::Long
  GetOptions($settings,qw(help tempdir=s downsample=s length=i sizeTo=s));
  die usage() if(@ARGV<1 || $$settings{help});
  # additional settings using ||= operator
  $$settings{tempdir}||=AKUtils::mktempdir(); # any temp file you make should go under the tempdir
  $$settings{downsample}||=0;
  die "ERROR: downsample should be a number and less than 1".usage() if($$settings{downsample}=~/[A-Za-z]/ || $$settings{downsample}>1);

  my $infile=[@ARGV];
  
  removeDuplicateReads($infile,$settings);

  return 0;
}

sub removeDuplicateReads{
  my($infile,$settings)=@_;
  my $linker="~" x 4; # tildes are probably not in the sequence, so it's safe for putting into the hash ID
  my $linkerFirstChar=substr($linker,0,1);
  my $l=$$settings{length}||0;

  for my $fastq(@$infile){
    my ($basename,$dir,$inSuffix) = fileparse($fastq,@IOextensions);

    # change the downsample parameter if settings{sizeTo} is set
    if($$settings{sizeTo}){
      $$settings{downsample}=findDownsamplingFromSizeto($fastq,$$settings{sizeTo},$settings);
    }

    my $poly=AKUtils::is_fastqPE($fastq,$settings)+1;
    if($inSuffix=~/\.fastq$/){
      open(FQ,'<',$fastq) or die "Could not open $fastq for reading: $!";
    }
    elsif($inSuffix=~/\.fastq\.gz/){
      open(FQ,"gunzip -c $fastq | ") or die "Could not open $fastq for reading: $!";
    }
    else{
      die "Could not determine the file type for reading based on your extension $inSuffix";
    }
    logmsg "Reading $fastq with poly=$poly";
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
        $hashId.="$linker$sequence2";
      }

      # if desired, only look at the first X nucleotides when considering a duplicate
      if($l){
        $hashId=~s/^(.{$l,$l}).*$linker(.+)/$1$linker$2/; # accept only X nucleotides from the front
        $hashId=~s/($linker.{$l,$l}).*($|$linker)/$1$2/g if($poly>1);
      }

      # downsampling at this step is more accurate than later, after dupes are removed
      if($$settings{downsample}){
        next if(rand(1)>$$settings{downsample});
      }
      
      $read{$hashId}=$read;
    }
    close FQ;

    logmsg "Writing unique reads from $fastq";
    while(my($hashId,$read)=each(%read)){
      print $read;
    }
  }
  logmsg "Done";
}

sub findDownsamplingFromSizeto{
  my($fastq,$sizeTo,$settings)=@_;

  # probably not the best way to do this, since run_assembly_readMetrics.pl could change col order
  my $origBp=`run_assembly_readMetrics.pl --fast '$fastq'|cut -f 3|tail -n 1`; chomp($origBp);
  my $downsample=$sizeTo/$origBp;
  $downsample=1 if($downsample > 1);
  logmsg "Estimating downsampling as $downsample of the original with requested size of $sizeTo";
  return $downsample;
}


sub usage{
  "Removes duplicate reads from a raw read set
   Usage: $0 read.fastq[.gz] > read.fastq
         $0 read.fastq[.gz] | gzip -c > read.fastq.gz
     --downsample 0.5  # filter 50% of reads out
     --length 100      # only consider up to 100bp when deciding if a read is a duplicate. Default: no limit
     -s downsample to the many bp. Internally, a new --downsample parameter is calculated
  ";
}

