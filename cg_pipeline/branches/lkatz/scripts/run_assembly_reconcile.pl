#!/usr/bin/env perl
# author: Lee Katz (lkatz@cdc.gov)

use strict;
use warnings;
use FindBin;
use lib "$FindBin::RealBin/../lib";
$ENV{PATH} = "$FindBin::RealBin:".$ENV{PATH};
use AKUtils;
use CGPipeline::Reconciliator;
use Getopt::Long;
use Data::Dumper;
use File::Basename;

$0=basename $0;
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help));
  die usage() if($$settings{help});
  
  my $reconciliator=CGPipeline::Reconciliator->new();
  
  for(@ARGV){
    $reconciliator->addVelvetAssembly({velvetDir=>$_});
  }
  print Dumper $reconciliator;
  
  return 0;
}

sub usage{
  "Reconciles two or more assemblies
  Usage: $0 *.assembly.fasta
  "
}
