#!/usr/bin/perl

# run_pipeline_NCBI: Turn your project into something submitable to NCBI
# Author: Lee Katz (lskatz@gatech.edu)

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
use Switch 'fallthrough';

my $verbose;

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH = $FSFind::LOG || *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

exit(main());
sub main(){
  $settings = AKUtils::loadConfig($settings);
  die usage() if(@ARGV<1);
  GetOptions($settings,qw(help! project=s tempdir=s keep!));
  die usage() if($$settings{help});
  $$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($$settings{keep}));

  # check for required files
  my @file=qw(assembly.fasta annotation.gb template.sbt);
  for(my $i=0;$i<@file;$i++){
    $file[$i]=File::Spec->rel2abs("$$settings{project}/$file[$i]");
    if(-s $file[$i] < 1){
      die "File does not exist or is 0 bytes: $file[$i]";
    }
  }
  my($assembly,$annotation,$template)=@file;

  my $template=createTemplate($settings);

  # add on anything necessary for the defline, etc
  $assembly=defineAssembly($assembly,$settings);
  logmsg "Tagged assembly file has been printed to $assembly";

  return 0;
}

# Creates an NCBI template file http://www.ncbi.nlm.nih.gov/genbank/tbl2asn2.html#Template
# Also see: http://www.ncbi.nlm.nih.gov/projects/HTGS/processing.html
# Uses data from config file
# Assumptions: first and last authors are the contacts
# TODO put in any necessary config options such as pipeline_authors_address
sub createTemplate{
  my($settings)=@_;
  
}

# Adds on additional tags to the assembly.
# Assumes a unique seqid is already present.
# Does not change the actual assembly data.
# TODO just grab the assembly from the genbank file to reduce the number of necessary files, to keep consistency, and to promote modularity/portability of CG-Pipeline
sub defineAssembly{
  my($assembly,$settings)=@_;
  my $outfile="$$settings{tempdir}/assembly.fasta";
  my @classification=split(/\s+/,$$settings{classification});
  my $scientificName=join(" ",@classification[0..1]);
  my $strain=$$settings{project};
  my $seqs=AKUtils::readMfa($assembly,$settings);
  my $newSeqs;
  while(my($seqname,$seq)=each(%$seqs)){
    # TODO define any other tags here.
    # Allowable tags: strain, isolate, chromosome, topology, locatino, molecule, technique, protein, genetic code
    if($seqname!~/organism=/){
      $seqname.=" [organism=$scientificName]";
    }
    if($seqname!~/strain=/){
      $seqname.=" [strain=$strain]";
    }
    $$newSeqs{$seqname}=$seq;
  }

  # make a new assembly file and return its new name
  AKUtils::printSeqsToFile($newSeqs,$outfile,$settings);
  return $outfile;
}

sub usage{
  "Turns your genome project directory into something that you can submit to NCBI
  Usage: perl $0 -p project
  Where the project directory has the following files:
    assembly.fasta
    annotation.gb
  Arguments
    -p project directory name
    -h (no arguments) this help menu
    -k (no arguments) to keep temporary directory around. 
      It will be kept without specifying if the temp dir name is given.
    -t temp directory name
      This is the directory where temporary files are stored. If specified,
      the temp directory will not be deleted.
  ";
}
