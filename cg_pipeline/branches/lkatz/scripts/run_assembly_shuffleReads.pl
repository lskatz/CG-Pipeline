#!/usr/bin/env perl
# shuffle or deshuffle sequences.  Good for fastq files only right now.
# Author: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(deshuffle help));
  die usage() if($$settings{help});

  my $numReads;
  if($$settings{deshuffle}){
    die "ERROR: need a file to deshuffle\n". usage() if(@ARGV<1);
    my($seqFile)=@ARGV;
    $numReads=deshuffleSeqs($seqFile,$settings);
  } else {
    die "ERROR: need two files to shuffle\n". usage() if(@ARGV<2);
    my($seqFile1,$seqFile2)=@ARGV;
    $numReads=shuffleSeqs($seqFile1,$seqFile2,$settings);
  }
  return 0;
}

sub deshuffleSeqs{
  my($seqFile,$settings)=@_;
  open(SHUFFLED,"<",$seqFile) or die "Could not open shuffled fastq file $seqFile: $!";
  my $i=0;
  while(<SHUFFLED>){
    my $mod=$i%8;
    print STDOUT $_ if($mod<4);
    print STDERR $_ if($mod>=4);
    
    $i++;
  }
  close SHUFFLED;
  my $numReads=$i/4;
  return $numReads;
}

sub shuffleSeqs{
  my($seqFile1,$seqFile2,$settings)=@_;
  open(MATE1,"<",$seqFile1) or die "Could not open $seqFile1:$!";
  open(MATE2,"<",$seqFile2) or die "Could not open $seqFile2:$!";
  my $i=0;
  while(my $out=<MATE1>){
    $out.=<MATE1> for(1..3);
    $out.=<MATE2> for(1..4);
    print STDOUT $out;
    $i++;
  }
  close MATE1; close MATE2;
  my $numReads=$i/4;
  return $numReads;
}

sub usage{
  local $0=fileparse($0);
  "Shuffle or deshuffle sequences
  Usage: $0 -d shuffled.fastq > file_1.fastq 2> file_2.fastq
  $0 file_1.fastq file_2.fastq > shuffled.fastq
  $0 file_[12].fastq > shuffled.fastq
    -d for deshuffle
  Due to the double redirection, error messages are hidden. A user should check the exit code to see if the program executed correctly.
  "
}
