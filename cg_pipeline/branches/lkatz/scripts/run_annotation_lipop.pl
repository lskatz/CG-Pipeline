#!/usr/bin/env perl

use strict;
use FindBin;
use lib "$FindBin::RealBin/../lib";
$ENV{PATH} = "$FindBin::RealBin:".$ENV{PATH};
use AKUtils qw(logmsg);
use CGPipelineUtils;
use File::Basename;
use Data::Dumper;
use Getopt::Long;

$0=fileparse $0;
exit(main());

sub main{
  my $settings={};
  $settings=AKUtils::loadConfig($settings);
  GetOptions($settings,qw(tempdir=s outfile=s));
  die usage() unless -f $ARGV[0];
  $$settings{tempdir}||=AKUtils::mktempdir();
  $$settings{outfile}||="$ARGV[0].lipop.sql";

  my $lipopOut="$$settings{tempdir}/lipop.out";
  my $lipop=AKUtils::fullPathToExec("LipoP");
  system("$lipop --nod -short < $ARGV[0] > $$settings{tempdir}/lipop.out");
  die if $?;

  open(LIPOP,"<",$lipopOut) or die "ERROR: could not open lipoP output file:$!";
  open(OUT,">",$$settings{outfile}) or die "ERROR: could not write to outfile: $!";
  print OUT "# SpI: signalpeptidaseI; SpII: lipoprotein signal peptide; TMH: n-terminal transmembrane helix; CYT: cytoplasmic\n# name location score margin CleavageStart CleavageEnd\n";
  while(<LIPOP>){
    s/^#//; # remove starting # which is there for no reason
    s/^\s*|\s*$//g; # trim whitespace
    my($name,$loc,$score,$margin,$cleavage)=split /\s+/;
    next if(!$margin);
    $score=~s/score=//;
    $margin=~s/margin=//;
    $cleavage=~s/cleavage=//;
    my($start,$end)=split(/-/,$cleavage);

    print OUT join("|",$name,$loc,$score,$margin,$start,$end)."\n";
  }
  close OUT;
  close LIPOP;
  return 0;
}

sub usage{
  local $0=fileparse $0;
  "Usage: $0 input.faa [-o outfile -t tempdir/]"
}
