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
use Bio::Perl;

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH = *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

exit(main());

sub main() {
  $settings = AKUtils::loadConfig($settings);

  die(usage($settings)) if @ARGV < 1;

  my @cmd_options = qw(keep tempdir=s outfile=s force assembly=s@ reads=s@ minContigSize=i);
  GetOptions($settings, @cmd_options) or die;

  $$settings{cpus}=AKUtils::getNumCPUs();
  $$settings{minContigSize}||=500;

  my @assembly=@{ $$settings{assembly} } or die "Need assemblies:\n".usage($settings);
  $$settings{outfile} ||= "$0.out.fasta"; # TODO determine if an ace should be the output
  $$settings{outfile} = File::Spec->rel2abs($$settings{outfile});
  open(FH, '>', $$settings{outfile}) or die("Error writing to output file $$settings{outfile}: $!");
  close FH;

  $$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($$settings{keep}));
  logmsg "Temporary directory is $$settings{tempdir}";

  if($$settings{reads}){
    my $reads=combineReadsFiles($settings);
    push(@assembly,$reads);
  }

  my $combinedAssembly=combineAllAssemblies(\@assembly,$settings);
  system("cp $combinedAssembly $$settings{outfile}"); 
  die "Could not copy final file $combinedAssembly to $$settings{outfile}" if $?;

  logmsg "Output file is " . $$settings{outfile};

  return 0;
}

sub combineReadsFiles{
  my($settings)=@_;
  $$settings{readsFile}="$$settings{tempdir}/extraReads.fasta";
  logmsg "Combining all reads into a single reads file at $$settings{readsFile}";
  my $reads={};
  for my $file(@{$$settings{reads}}){
    my $seqs=AKUtils::readMfa($file,$settings);
    for my $s(keys(%$seqs)){
      $$reads{$s}=$$seqs{$s};
    }
  }
  AKUtils::printSeqsToFile($reads,$$settings{readsFile},$settings);
  return $$settings{readsFile};
}

sub combineAllAssemblies{
  my($assembly,$settings)=@_;
  
  logmsg "Sorting the assemblies by combining assembly metrics";
  my @assembly=sort({
    assemblyScore($b) <=> assemblyScore($a);
  } @$assembly);
  $assembly=\@assembly;

  my $combinedAssembly=shift(@$assembly); # the first assembly is the first "combined" assembly
  for(my $i=0;$i<@$assembly;$i++){
    $$settings{is_reads_file}=1 if($i+1==@$assembly && $$settings{readsFile});
    my $otherAssembly=$$assembly[$i];
    next if(-s $otherAssembly < 1);
    $combinedAssembly=combine2Assemblies($combinedAssembly,$otherAssembly,$settings);
  }

  logmsg "Filtering sequences by contig size of $$settings{minContigSize}";
  my $filteredSeqs;
  my $seqs=AKUtils::readMfa($combinedAssembly,$settings);
  while(my($id,$sequence)=each(%$seqs)){
    $$filteredSeqs{$id}=$sequence if(length($sequence)>$$settings{minContigSize});
  }
  $$settings{order_seqs_by_length}=1;
  AKUtils::printSeqsToFile($filteredSeqs,$combinedAssembly,$settings);
  
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

  # combine the input fastas into one, with unique IDs
  #system("cat '$a1' '$a2' > $combined_fasta_file");
  #die "Problem with creating input file" if $?;
  my @seq;
  for my $file($refGenome,$queryGenome){
    my $in=Bio::SeqIO->new(-file=>$file);
    while(my $seq=$in->next_seq){
      $seq->id(join("_",$file,$seq->id));
      push(@seq,$seq);
    }
  }
  my $seqout=Bio::SeqIO->new(-file=>">$combined_fasta_file");
  $seqout->write_seq(@seq);

  logmsg "Running Minimus2 with reference genome $refGenome and query genome $queryGenome";
  system("toAmos -s '$combined_fasta_file' -o '$$settings{tempdir}/minimus.combined.afg'");
  die "Problem with toAmos with command\n  toAmos -s '$combined_fasta_file' -o '$$settings{tempdir}/minimus.combined.afg'" if $?;
  system("minimus2 -D REFCOUNT=$numContigs '$$settings{tempdir}/minimus.combined'");
  die "Problem with Minimus2" if $?;

  # recover singletons that pass the filter
  my %allseqs=(); my $i=0;
  for my $file("$$settings{tempdir}/minimus.combined.fasta","$$settings{tempdir}/minimus.combined.singletons.seq"){
    my $seqs=AKUtils::readMfa($file,$settings);
    while(my($id,$sequence)=each(%$seqs)){
      $allseqs{++$i}=$sequence if (length($sequence)>=$$settings{minContigSize});
    }
  }
  AKUtils::printSeqsToFile(\%allseqs,"$$settings{tempdir}/contigsAndSingletons.fasta");

  return "$$settings{tempdir}/contigsAndSingletons.fasta";
}

# make a crude assembly score for sorting the best assemblies
sub assemblyScore{
  my($a,$settings)=@_;
  my $tmp= `run_assembly_metrics.pl $a -s assemblyScore -n`;
  chomp($tmp);
  return $tmp;
}

sub assembly_metrics{
  my($a,$settings)=@_;
  my %seqMetric;
  return %seqMetric if(-s $a < 1);
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
    -a assembly file in fasta format
      multiple assemblies are allowed (and expected)
    -r reads file in fasta format
      Reads that have not been used in the individual assemblies
    -o outputfile
      assembly output fasta file
    optional parameters:
    -m minimum size of a contig
    -t tempdir
    -f
      force
    -k
      keep temp files
  "
}
