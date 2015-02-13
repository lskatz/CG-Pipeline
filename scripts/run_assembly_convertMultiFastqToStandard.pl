#!/usr/bin/env perl

# Some fastq files are nonstandard, and each entry spans more than four rows, making them difficult to parse.
# This script converts those entries.
# Author: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename;

sub logmsg{local $0=fileparse $0; print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(infile=s outfile=s minLength=i help));
  die usage() if($$settings{help});
  my $in=$$settings{infile} or logmsg "No infile. Reading from stdin";
  my $out=$$settings{outfile} or logmsg "No outfile. Printing to stdout";

  convert($in,$out,$settings);

  return 0;
}

sub convert{
  my($in,$out,$settings)=@_;

  # Where to get input from?
  my $infp;
  if($in){
    open($infp,$in) or die "Could not open $in: $!";
  } else {
    $infp=*STDIN;
  }

  # Where to print to?
  my $outfp;
  if($out){
    open($outfp,">$out") or die "Could not open $out for writing: $!";
  } else {
    $outfp=*STDOUT;
  }

  while(my $id=<$infp>){
    my($qual,$sequence);
    chomp($id);
    $id=~s/^@//;
    while(my $tmp=<$infp>){
      last if($tmp=~/^\+(\Q$id\E)?/);
      $sequence.=$tmp;
    }
    my $seqLength=length($sequence);

    read($infp,$qual,$seqLength) or die "Problem reading $seqLength characters from $in because $!";
    die "For sequence $id, I got a qual that is a different length than intended\n$qual" if(length($qual)!=$seqLength);
    #for(1..$seqLength){
    #  $qual.=getc(IN);
    #}
    $qual=~s/\s+//g;
    $sequence=~s/\s+//g;

    next if($$settings{minLength} && length($sequence)<$$settings{minLength});

    print $outfp "\@$id\n$sequence\n+\n$qual\n";
  }
  close $outfp;
  close $infp;

  return 1;
}

sub usage{
  "Converts a multi-line fastq file to the standard 4-lines per entry, standard form
  Usage: $0 -i nonstandard.fastq -o standard.fastq [other options]
  -m minimum length for a sequence to be outputted, in base pairs
  ";
}
