#!/usr/bin/env perl
# Downsample a fastq file

use strict;
use warnings;
use Getopt::Long;
use FindBin;
use File::Basename qw/basename/;

use lib "$FindBin::RealBin/../lib";
use AKUtils qw/logmsg/;

$AKUtils::LOG=*STDERR;

$0=basename($0);

exit(main());
sub main{
  my $settings={};
  GetOptions($settings,qw(downsample=f help)) or die $!;
  die usage() if($$settings{help});
  $$settings{downsample}//=1;
  die "ERROR: downsampling must be between 0 and 1" if($$settings{downsample} > 1 || $$settings{downsample} < 0);
  my $fastq=$ARGV[0];
  die usage() if(!$fastq);

  my $numReads=0;
  my $numReadsWritten=0;
  open(FASTQ,"gunzip -c $fastq |") or die "ERROR: could not open $fastq for reading: $!";
  while(my $peRead=<FASTQ>){
    $peRead.=<FASTQ> for(2..8); # read the remaining lines
    $numReads++;

    # Only accept the read if the downsample param is higher than a random decimal
    next if(rand() > $$settings{downsample});

    # Write the read
    print $peRead;
    $numReadsWritten++;

    # Progress report
    if($numReadsWritten % 100000 ==0){
      my $percentWritten=int($numReadsWritten/$numReads * 100);
      logmsg "$numReadsWritten reads written out of $numReads so far ($percentWritten%)";
    }
  }
  close FASTQ;

  # Final report
  my $percentWritten=int($numReadsWritten/$numReads * 100);
  logmsg "Finished writing $numReadsWritten out of $numReads ($percentWritten%)";

  return 0;
}

sub usage{
  "$0: Downsamples a paired end fastq.gz file
  Usage: $0 --downsample 0.1 in.fastq.gz > out.fastq
  "
}
