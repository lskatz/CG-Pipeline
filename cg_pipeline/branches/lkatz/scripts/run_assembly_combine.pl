#!/usr/bin/env perl

# run-assembly-reconciliation: combine assemblies
# Author: Lee Katz <lskatz@gatech.edu>

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
use File::Temp ('tempdir');
use File::Path;
use File::Spec;
use File::Copy;
use File::Basename;
use List::Util qw(min max sum shuffle);
use CGPipelineUtils;
use Data::Dumper;

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH = *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

exit(main());

sub main() {
  $settings = AKUtils::loadConfig($settings);

  die(usage()) if @ARGV < 1;

  my @cmd_options = qw(keep tempdir=s outfile=s force assembly=s@);
  GetOptions($settings, @cmd_options) or die;

  $$settings{cpus}=AKUtils::getNumCPUs();

  my @assembly=@{ $$settings{assembly} } or die "Need assemblies:\n".usage();
  $$settings{outfile} ||= "$0.out.fasta"; # TODO determine if an ace should be the output
  $$settings{outfile} = File::Spec->rel2abs($$settings{outfile});
  open(FH, '>', $$settings{outfile}) or die("Error writing to output file $$settings{outfile}: $!");
  close FH;

  $$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($$settings{keep}));
  logmsg "Temporary directory is $$settings{tempdir}";

  my $combinedAssembly=combineAllAssemblies(\@assembly,$settings);
  system("cp $combinedAssembly $$settings{outfile}"); 
  die "Could not copy final file $combinedAssembly to $$settings{outfile}" if $?;

  logmsg "Output file is " . $$settings{outfile};

  return 0;
}

sub combineAllAssemblies{
  my($assembly,$settings)=@_;
  my $combinedAssembly=shift(@$assembly);
  for(my $i=0;$i<@$assembly;$i++){
    my $otherAssembly=$$assembly[$i];
    $combinedAssembly=combine2Assemblies($combinedAssembly,$otherAssembly,$settings);
  }
  return $combinedAssembly;
}

sub combine2Assemblies{
  my($a1,$a2,$settings)=@_;
  my %metrics1=assembly_metrics($a1,$settings);
  my %metrics2=assembly_metrics($a2,$settings);
  my($numContigs1,$numContigs2)=($metrics1{numContigs},$metrics2{numContigs});
  my $combined_fasta_file = "$$settings{tempdir}/combined_in.fasta";

  my $numContigs=$numContigs1;
  my $refGenome=$a1;
  my $queryGenome=$a2;
  system("cat '$a1' '$a2' > $combined_fasta_file");
  die "Problem with creating input file" if $?;
  if($numContigs2>$numContigs1){
    $numContigs=$numContigs2;
    $refGenome=$a2;
    $queryGenome=$a1;
    system("cat '$a2' '$a1' > $combined_fasta_file");
    die "Problem with creating input file" if $?;
  }
  logmsg "Running Minimus2 with reference genome $refGenome and query genome $queryGenome";
  system("toAmos -s '$combined_fasta_file' -o '$$settings{tempdir}/minimus.combined.afg'");
  die "Problem with toAmos with command\n  toAmos -s '$combined_fasta_file' -o '$$settings{tempdir}/minimus.combined.afg'" if $?;
  system("minimus2 -D REFCOUNT=$numContigs '$$settings{tempdir}/minimus.combined'");
  die "Problem with Minimus2" if $?;
  return "$$settings{tempdir}/minimus.combined.fasta";
}

sub assembly_metrics{
  my($a,$settings)=@_;
  my %seqMetric;
  my $metrics=`run_assembly_metrics.pl '$a'`;
  for my $line(split(/\n/,$metrics)){
    my($key,$value)=split(/\t/,$line);
    $seqMetric{$key}=$value;
  }
  return %seqMetric;
}

sub usage{
  "Combines two or more assemblies into one assembly, using Minimus2
  Usage: $0 -a assembly1 -a assembly2 [-a ...] -o output.assembly.fasta
    -a assembly file
      at least two assembly files are needed.  You can also decide to use a reads file at your own risk.
    -o outputfile
      assembly output fasta file
    optional parameters:
    -t tempdir
    -f
      force
    -k
      keep temp files

  "
}
