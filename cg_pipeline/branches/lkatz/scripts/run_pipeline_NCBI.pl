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
use Bio::Perl;

my $verbose;

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH = $FSFind::LOG || *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

exit(main());
sub main(){
  $settings = AKUtils::loadConfig($settings);
  die usage() if(@ARGV<1);
  GetOptions($settings,qw(help! genbank=s tempdir=s keep!));
  die usage() if($$settings{help});
  $$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($$settings{keep}));

  # check for required files
  my $genbank=$$settings{genbank} or die "Error: missing genbank parameter\n".usage();
  $genbank=File::Spec->rel2abs($genbank);

  # add on anything necessary for the defline, etc
  my $assembly=defineAssembly($genbank,$settings);
  logmsg "Tagged assembly file has been printed to $assembly";

  # Need to figure out how to create ASN file
  # http://biostar.stackexchange.com/questions/10748/read-writer-for-asn
  #my $template=createTemplate($settings);

  my $tbl=createTbl($genbank,$settings);

  return 0;
}


# See: http://www.ncbi.nlm.nih.gov/Sequin/table.html#Table Layout
sub createTbl{
  my($genbank,$settings)=@_;
  die "TODO: this subroutine";
}

# Creates an NCBI template file http://www.ncbi.nlm.nih.gov/genbank/tbl2asn2.html#Template
# Also see: http://www.ncbi.nlm.nih.gov/projects/HTGS/processing.html
# Uses data from config file
# Assumptions: first and last authors are the contacts
# TODO put in any necessary config options such as pipeline_authors_address
sub createTemplate{
  my($settings)=@_;

  # make array of authors
  my @author;
  my $a=1;
  while($$settings{"author$a"}){
    my %a;
    ($a{name},$a{affiliation},$a{div},$a{city},$a{country},$a{street},$a{email},$a{phone},$a{zipCode})=split(/\|/,$$settings{"author$a"});
    push(@author,\%a);
    $a++;
  };
  print Dumper \@author;exit;

  # make corresponding authors

  # make comprehensive author list
}


# Adds on additional tags to the assembly.
# Assumes a unique seqid is already present.
# Does not change the actual assembly data.
# TODO just grab the assembly from the genbank file to reduce the number of necessary files, to keep consistency, and to promote modularity/portability of CG-Pipeline
sub defineAssembly{
  my($genbank,$settings)=@_;
  my $outfile="$$settings{tempdir}/assembly.fasta";
  my $seqin=Bio::SeqIO->new(-file=>$genbank,-format=>"genbank");
  my $seqout=Bio::SeqIO->new(-file=>">$outfile");

  while(my $seq=$seqin->next_seq){
    #TODO use Bio::Taxon instead because Bio::Species is becoming deprecated

    #TODO there is not a reliable way of getting the genome name yet
    # Get the strain name from the locus name, e.g. 2010V-1031_0001 => 2010V-1031
    my $scientificName=$seq->species->binomial;
    my @f=split(/_/,$seq->id);
    my $seqid=pop(@f);
    my $strain=join("_",@f);

    $seq->desc("[organism=$scientificName] [strain=$strain]");
    $seq->id($seqid);

    $seqout->write_seq($seq);

    # TODO get classification and scientific name from the genbank file
    # TODO define any other tags here.
    # Allowable tags: strain, isolate, chromosome, topology, locatino, molecule, technique, protein, genetic code
    #my @feat=$seq->get_SeqFeatures;
  }

  # make a new assembly file and return its new name
  return $outfile;
}

sub usage{
  "Turns your genome project directory into something that you can submit to NCBI
  Usage: perl $0 -g genbank.gbk
  Arguments
    -g genbank annotated file
    -h (no arguments) this help menu
    -k (no arguments) to keep temporary directory around. 
      It will be kept without specifying if the temp dir name is given.
    -t temp directory name
      This is the directory where temporary files are stored. If specified,
      the temp directory will not be deleted.
  ";
}
