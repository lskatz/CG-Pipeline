#!/usr/bin/env perl

# run-assembly-convertOpgen: convert opgen data into paired end reads
# Author: Lee Katz (Lkatz@cdc.gov)

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
use XML::LibXML::Reader;
use Bio::Perl;

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH = $FSFind::LOG || *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

exit(main());

sub main() {
	$settings = AKUtils::loadConfig($settings);
	die usage() if(@ARGV<1);

	my @cmd_options = qw(infile=s keep! tempdir=s outfile=s restrictionSequence=s);
	GetOptions($settings, @cmd_options) or die;
  
  my $restrictionMapDocument=$$settings{infile} or die "Error: missing infile\n".usage();
  my $outfile=$$settings{outfile} || "$0.fastq";
  my $restrictionSequence=$$settings{restrictionSequence} or die "Error: missing restriction enzyme sequence\n".usage();
  $$settings{tempdir} ||= tempdir($$settings{tempdir} or File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($$settings{keep}));
  $$settings{qualScore}||=60;

  my @restrictionFragmentLength=readMap($restrictionMapDocument,$settings);
  logmsg "Found ".scalar(@restrictionFragmentLength)." sites";
  
  # make a raw reads file
  my $reads=writeReads($restrictionSequence,\@restrictionFragmentLength,$settings);

  move($reads,$outfile);
  logmsg "reads are in $outfile";

  return 0;
}

sub writeReads{
  my($restrictionSequence,$restrictionFragmentLength,$settings)=@_;
  #die "TODO: make pair of files with the distance between the reads encoded in the filename. Consequently, make a directory of files";
  
  my $dir=$$settings{tempdir};
  my %outseq;
  my $restrictionQual="$$settings{qualScore} " x length($restrictionSequence);
  chomp($restrictionQual);

  my $i=0;
  for my $rLength(@$restrictionFragmentLength){
    if(!$outseq{'f'.$rLength}){
      $outseq{'f'.$rLength}=Bio::SeqIO->new(-format=>"fastq-illumina",-file=>">$dir/f$rLength.fastq");
      $outseq{'r'.$rLength}=Bio::SeqIO->new(-format=>"fastq-illumina",-file=>">$dir/r$rLength.fastq");
    }
    #die "TODO make 20 reads per site to ensure confidence during the assembly";
    my $readid="restrictionsite".$i++;
    my $seq=Bio::Seq::Quality->new(-id=>$readid,-seq=>$restrictionSequence,-qual=>$restrictionQual);
    $outseq{'f'.$rLength}->write_seq($seq);
    $outseq{'r'.$rLength}->write_seq($seq);
  }
  return $dir;
}

sub readMap{
  my($map,$settings)=@_;
 
  my @rLength;

  my $reader = new XML::LibXML::Reader(location => $map) or die "cannot read $map\n";
  while($reader->read){
    if($reader->name eq 'RESTRICTION_MAP' && $reader->getAttribute('INSILICO') eq 'false'){
      while($reader->read){
        last if($reader->name eq 'RESTRICTION_MAP'); # break out at the closing tag
        next if($reader->name ne 'F'); # F for fragment
        push(@rLength,$reader->getAttribute('S'));
      }
    }
  }

  return @rLength;
}

sub processNode {
      my $reader = shift;
      printf "%d %d %s %d\n", ($reader->depth,
                               $reader->nodeType,
                               $reader->name,
                               $reader->isEmptyElement);
  }


sub usage{
  "Create a paired end read input file for assembly, from OpGen data
  Usage: $0 -i opgen.xml -o outfile -r sequence
    -i input file
      This is XML as given by OpGen or MapSolver
    -o outfile
      This is your paired end read output file
    -t temporary directory
      This is where temporary files are stored
    -r sequence
      Restriction enzyme sequence (e.g. GCTAGC for NheI)
    -h 
      this help menu
  ";
}

