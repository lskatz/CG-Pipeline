#!/usr/bin/env perl

# Predict CRISPRs in a fasta file
# Author: Lee Katz <lkatz@cdc.gov>
# TODO output the entire crispr in a GFF line and not just DRs

use strict;
use warnings;
use FindBin;
use lib "$FindBin::RealBin/../lib";
$ENV{PATH} = "$FindBin::RealBin:".$ENV{PATH};
use AKUtils qw/logmsg/;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use File::Copy qw/copy/;
use Bio::AlignIO;

local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
exit(main());

sub main{
  local $0=fileparse $0;
  my $settings={
    appname=>"cgpipeline",
    minCrisprSize=>23,
    maxCrisprSize=>55,
  };
  GetOptions($settings,qw(outfile=s help windowSize=i stepSize=i tempdir=s minCrisprSize=i maxCrisprSize=i verbose));
  $$settings{tempdir}||=AKUtils::mktempdir($settings);

  $$settings{pilerCr}=AKUtils::fullPathToExec("pilercr");
  my $outfile=$$settings{outfile} || "$0.gff";
  die usage() if(@ARGV<1);
  my @assembly=@ARGV;
  
  my @gff;
  for my $assembly(@assembly){
    my $pilerCrOut=runPilerCr($assembly,$settings);
    my $gff=pilerCrToGff($pilerCrOut,$settings);
    push(@gff,$gff);
  }
  
  system("cat ".join(" ",@gff)." > $outfile");
  die if $?;

  return 0;
}

sub runPilerCr{
  my($assembly,$settings)=@_;

  logmsg "$assembly";
  my $basename=fileparse $assembly;
  my $pilerCrOut="$$settings{tempdir}/$basename.pilerCr.out";
  my $pilerCrDr="$$settings{tempdir}/$basename.pilerCr.dr.fasta";
  # pilercr -in P_aeruginoasP14.fasta -out /dev/stdout -seq pilecr.seq -noinfo
  system("$$settings{pilerCr} -in '$assembly' -out '$pilerCrOut' -seq '$pilerCrDr' -noinfo -quiet");
  die if $?;

  return ($pilerCrOut,$pilerCrDr) if wantarray;
  return $pilerCrOut;
}

sub pilerCrToGff{
  my($pilerCrOut,$settings)=@_;
  my $gff="$pilerCrOut.gff";
  open(PILERCR,"<",$pilerCrOut) or die "Could not open piler-cr file $pilerCrOut:$!";
  open(GFF,">",$gff) or die "Could not open $gff for writing:$!";
  my $stage="";
  while(my $line=<PILERCR>){
    if($line=~/DETAIL REPORT/){
      $stage="detail";
    }
    elsif($line=~/SUMMARY BY SIMILARITY/){
      $stage="similarity";
    }
    elsif($line=~/SUMMARY BY POSITION/){
       $stage="position"; 
    }
    # detailed DRs
    elsif($stage eq "detail" && $line=~/Array\s+(\d+)/){
      my $crisprId=$1;
      my $seqname=<PILERCR>; chomp($seqname);
      $seqname=~s/\s+.*//; # remove anything after whitespace
      my $undef=<PILERCR> for(1..3); # discard: blank, headers, ===
      while(<PILERCR>){
        last if /===/;
        my @F=split /\s+/; @F=grep(!/^\s*$/,@F);
        print GFF join("\t",$seqname,"Piler-CR","direct_repeat",$F[0],($F[0]+$F[1]-1),$F[2],'.','.',"")."\n";
      }
    }

  }
  close GFF;
  close PILERCR;
  return $gff;
}

sub usage{
  my($settings)=@_;
  local $0=fileparse $0;
  "Finds CRISPRs in a fasta file
  Usage: $0 assembly.fasta -o outfile.gff
  "
}
