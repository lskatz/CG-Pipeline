#!/usr/bin/env perl

use strict;
use warnings;

exit(main(@ARGV));

sub main{
  my($sFastq, $prefix)=@_;
  die usage() if(!defined($prefix));

  my($mate1,$mate2)=($prefix."_1.fastq",$prefix."_2.fastq");
  open(SHUFFLED,$sFastq) or die "Could not open shuffled fastq file $sFastq: $!";
  open(MATE1,">",$mate1) or die "Could not open mate1 for writing: $!";
  open(MATE2,">",$mate2) or die "Could not open mate1 for writing: $!";
  my $i=0;
  while(<SHUFFLED>){
    my $mod=$i%8;
    print MATE1 $_ if($mod<4);
    print MATE2 $_ if($mod>=4);
    
    $i++;
  }
  close SHUFFLED; close MATE1; close MATE2;
  return 0;
}

sub usage{
  "Deshuffle a fastq file
  Usage: $0 shuffled.fastq prefix
    Output files will be prefix_1.fastq and prefix_2.fastq
  "
}
