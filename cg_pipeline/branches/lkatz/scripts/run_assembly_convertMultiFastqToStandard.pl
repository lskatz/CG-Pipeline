#!/usr/bin/env perl

# Some fastq files are nonstandard, and each entry spans more than four rows, making them difficult to parse.
# This script converts those entries.
# Author: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(infile=s outfile=s minLength=i));
  my $in=$$settings{infile} or die usage();
  my $out=$$settings{outfile} or die usage();

  convert($in,$out,$settings);

  return 0;
}

sub convert{
  my($in,$out,$settings)=@_;
  open(IN,$in) or die "Could not open $in: $!";
  open(OUT,">$out") or die "Could not open $out for writing: $!";
  while(my $id=<IN>){
    my($qual,$sequence);
    chomp($id);
    $id=~s/^@//;
    while(my $tmp=<IN>){
      last if($tmp=~/^\+(\Q$id\E)?/);
      $sequence.=$tmp;
    }
    my $seqLength=length($sequence);

    for(1..$seqLength){
      $qual.=getc(IN);
    }
    $qual=~s/\s+//g;
    $sequence=~s/\s+//g;

    next if($$settings{minLength} && length($sequence)<$$settings{minLength});

    print OUT "\@$id\n$sequence\n+\n$qual\n";
  }
  close OUT;
  close IN;

  return 1;
}

sub usage{
  "Converts a multi-line fastq file to the standard 4-lines per entry, standard form
  Usage: $0 -i nonstandard.fastq -o standard.fastq [other options]
  -m minimum length for a sequence to be outputted, in base pairs
  ";
}
