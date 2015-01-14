#!/usr/bin/env perl

# run-prediction-prodigal: Perform gene prediction with prodigal on a fasta file
# Author: Lee Katz (lkatz@cdc.gov)

package PipelineRunner;
my ($VERSION) = ('$Id: $' =~ /,v\s+(\d+\S+)/o);

my $settings = {
  appname => 'cgpipeline',
};
my $stats;

use strict;
use FindBin;
use lib "$FindBin::RealBin/../lib";
$ENV{PATH} = "$FindBin::RealBin:".$ENV{PATH};
use AKUtils qw(logmsg);

use Getopt::Long;
use File::Path;
use File::Spec;
use File::Copy;
use File::Basename;
use Bio::Seq;
use Bio::SeqIO;
#use Bio::SeqFeature::Gene::GeneStructure;
use AKUtils;
use Data::Dumper;

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH = *STDERR; print $FH "$0: ".(caller(1))[3].": @_\n";}

exit(main());

sub main() {
  $settings = AKUtils::loadConfig($settings);
  GetOptions($settings,qw(help tempdir=s keep));
  die usage() if(!@ARGV || $$settings{help});
  $$settings{tempdir}||=AKUtils::mktempdir($settings);

  # check for each assembly
  my @assembly=@ARGV;
  for(@assembly){
    die "ERROR: could not locate assembly file $_" if(!-f $_);
    #$_=File::Spec->rel2abs($_);
  }
  # concat all assemblies into one
  my $assembly="$$settings{tempdir}/assembly.fasta";
  system("cat ".join(" ",@assembly)." > $assembly"); die if $?;

  # find where prodigal is installed
  $$settings{exec}=AKUtils::fullPathToExec("prodigal");

  logmsg "Training Prodigal";
  my $training=train($assembly,$settings);
  logmsg "Predicting with Prodigal";
  my $gff=runPrediction($assembly,$training,$settings);

  system("cat '$gff'"); die if $?;

  return 0;
}

# get the training data from prodigal
sub train{
  my($assembly,$settings)=@_;
  my $training="$$settings{tempdir}/prodigal.training";
  return $training if(-f $training);
  my $prodigal=$$settings{exec};
  #prodigal -i assembly.fasta -t prodigal.training 
  system("$prodigal -i $assembly -t $training"); die if $?;
  return $training;
}

sub runPrediction{
  #prodigal -i assembly.fasta -t prodigal.training -f gff > prodigal.gff
  my($assembly,$training,$settings)=@_;
  my $prodigal=$$settings{exec};
  my $gff="$$settings{tempdir}/prodigal.gff";
  system("$prodigal -i $assembly -t $training -f gff > $gff"); die if $?;
  return $gff; 
}

sub usage{
  "Predicts for coding sequences in an assembly using Prodigal
  Usage: $0 assembly.fasta [assembly2.fasta] > prediction.gff
    -t tempdir/
  "
}
