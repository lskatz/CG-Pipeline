#!/usr/bin/env perl

# run-assembly-opgen: incorporate an opgen xml to make an assembly better
# Author: Lee Katz (Lkatz@cdc.gov)

package PipelineRunner;
my ($VERSION) = ('$Id: $' =~ /,v\s+(\d+\S+)/o);

my $settings = {
	appname => 'cgpipeline',
};
$|++;

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
use Bio::Assembly::IO;
use Bio::Assembly::Scaffold;
use Math::Round qw(round nearest);

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH = *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

exit(main());

sub main() {
	$settings = AKUtils::loadConfig($settings);
	die usage() if(@ARGV<1);

	my @cmd_options = qw(restrictionXml=s ace=s sequenceReads=s qualityReads=s keep! help! tempdir=s outfile=s);
	GetOptions($settings, @cmd_options) or die;
  die usage() if($$settings{help});
  
  my $restrictionMapDocument=$$settings{restrictionXml} or die "Error: missing xml infile\n".usage();
  my $outfile=$$settings{outfile} || "$0.fasta";
  my $ace=$$settings{ace} or die "Error: missing input ace file\n".usage();
  my $reads=$$settings{sequenceReads} or die "Error: missing sequence reads\n".usage();
  my $qual=$$settings{qualityReads} or die "Error: missing quality reads\n".usage();
  $$settings{tempdir} ||= tempdir($$settings{tempdir} or File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($$settings{keep}));

  # find which contigs are not mapped or are mapped more than once
  my $wellMappedContigs=categorizeContigs($restrictionMapDocument,$settings);

  # make a scaffold of mapped contigs
  my ($referenceGenome,$placedReads)=makeAssemblyHack($wellMappedContigs,$ace,$settings);
  logmsg "Created a reference genome of well-placed contigs $referenceGenome and noted which reads were assembled";

  # make fastq input of reads not found in the reference genome
  my $unplacedReads=makeReadsInput($reads,$qual,$placedReads,$settings);
  logmsg "Created a file with unplaced reads $unplacedReads";

  # run reference assembly of the unplaced reads vs the reference genome
  logmsg "Running mapping assembly";
  my $newbler_run=runNewblerMapping([$unplacedReads],[$referenceGenome],$settings);
  my $combined=combineNewblerMapping([$referenceGenome,"$newbler_run/mapping/454AllContigs.fna"],$settings);

  # TODO check the reference assembly output to see if there is still a misassembly. This might have to be a manual step for now.

  move($combined,$outfile);
  logmsg "Assembly is in $outfile";

  return 0;
}

# the reference genome should be the first one given
sub combineNewblerMapping{
  my($assemblies,$settings)=@_;
  $$settings{minCoverage}||=10;

  # concatenate and rename IDs to ensure unique seqids
  my $idCounter;
  my $minimusInfile="$$settings{tempdir}/minimus.combined.fna";
  my $seqout=Bio::SeqIO->new(-file=>">$minimusInfile");
  for my $a(@$assemblies){
    my $seqin=Bio::SeqIO->new(-file=>$a);
    while(my $seq=$seqin->next_seq){
      # don't accept a super-small coverage depth
      if($seq->desc =~ /numreads=(\d+)/){
        next if($1 < $$settings{minCoverage});
      }
      $seq->id(++$idCounter);
      $seqout->write_seq($seq);
    }
  }

  my $numContigs=`grep -c ">" $$assemblies[0]` + 0;
  system("toAmos -s '$minimusInfile' -o '$$settings{tempdir}/minimus.combined.afg'"); die if $?;
  system("minimus2 -D REFCOUNT=$numContigs '$$settings{tempdir}/minimus.combined'"); die if $?;
  return "$$settings{tempdir}/minimus.combined.fasta";
}

sub runNewblerMapping($$$) {
  my ($input_files, $ref_files, $settings) = @_;
  my $run_name = "$$settings{tempdir}/P__runMapping";
  #return $run_name; #debugging;
  logmsg "Executing Newbler mapping project $run_name";

  system("newMapping '$run_name'"); die if $?;
  foreach my $ref_file (@$ref_files) {
    system("setRef '$run_name' '$ref_file'"); die if $?;
  }
  foreach my $sff_file (@$input_files) {
    system("addRun '$run_name' '$sff_file'"); die if $?;
  }
  system("runProject '$run_name'"); die if $?;
  return $run_name;
}

sub makeReadsInput{
  my($reads,$qual,$placedReads,$settings)=@_;
  my (%placedReads,@unplacedReads);
  my $readsFile="$$settings{tempdir}/unplacedReads.fasta";
  $placedReads{$_}++ for (@$placedReads); # index the placed reads

  # make quality sequence hash
  logmsg "Creating the hash of quality sequences";
  my $qualin=Bio::SeqIO->new(-file=>$qual);
  my %qual;
  while(my $seq=$qualin->next_seq){
    $qual{$seq->id}=$seq;
  }

  logmsg "Writing reads file and quality file for unassembled reads";
  my $readsIn=Bio::SeqIO->new(-file=>$reads);
  my $readsOut=Bio::SeqIO->new(-file=>">$readsFile");
  my $qualOut=Bio::SeqIO->new(-file=>">$readsFile.qual");
  while(my $seq=$readsIn->next_seq){
    next if($placedReads{$seq->id});
    $readsOut->write_seq($seq);
    $qualOut->write_seq($qual{$seq->id});
  }
  return $readsFile;
}

sub makeAssemblyHack{
  my($placedContig,$ace,$settings)=@_;
  my(@contig,@line,@sequence,@placedReads);

  my $i=0; #contig counter
  my $lineNumber=0;
  open(ACE,$ace) or die "Could not open $ace because $!\n";
  CONTIG:
  while(<ACE>){
    if(/^CO/){
      $contig[$i]=$_;
      $contig[$i]=~s/^CO\s+(\S+).*$/$1/;
      $contig[$i]=~s/^\s+|\s+$//g;
      if(!$$placedContig{$contig[$i]}){
        delete($contig[$i]);
        # burn through lines until the contig is over
        while(<ACE>){
          if(/^AF/){ #burn through reads too
            while(<ACE>){
              next CONTIG if(/^\s*$/);
            }
          }
        }
        next;
      }
      while(<ACE>){
        chomp;
        s/\*+//g;
        $sequence[$i].=$_;
        if(/^\s*$/){
          $sequence[$i]=~s/(.{60})/$1\n/g;
          last;
        }
      }
      $i++; # next contig
    }
    elsif(/^AF/){
      chomp;
      s/^AF\s+(\S+).*$/$1/;
      push(@placedReads,$_);
    }
  }

  close ACE;

  my $assemblyOut="$$settings{tempdir}/wellPlacedContigs.fasta";
  open(OUT,">$assemblyOut") or die "Could not open the temporary assembly file $assemblyOut because $!\n";
  for(my $i=0;$i<@contig;$i++){
    print OUT ">$contig[$i]\n$sequence[$i]\n";
  }
  close OUT;
  return ($assemblyOut,\@placedReads);
}

sub makeAssembly{
  my($placedContig,$ace,$settings)=@_;

  my @keptContig;
  my $aIO=Bio::Assembly::IO->new(-file=>$ace);
  print Dumper $aIO;
  my $assembly=$aIO->next_assembly;
  exit;
  print Dumper $assembly;exit;
  foreach my $contig($assembly->all_contigs){
    print $contig->id." <=\n";
    if($$placedContig{$contig->id}){
      print Dumper $contig;exit;
      push(@keptContig,$contig);
    }
  }
  print Dumper \@keptContig;exit;
  my $assemblyout=Bio::Assembly::Scaffold->new();
}

sub categorizeContigs{
  my($map,$settings)=@_;
  
  # Count contig names.
  # Each contig will get counted if it appears on the XML anywhere.
  # Each contig will get counted if it maps to the in vitro genome.
  # Contigs that are mapped twice will have a count > 2.
  # Contigs that are mapped once will have a count < 2.
  # Well-mapped contigs will have a count == 2.
  my %contig;
  my $reader = new XML::LibXML::Reader(location => $map) or die "cannot read $map\n";
  while($reader->read){
    if($reader->name =~ /^(MAP_ALIGNMENT|RESTRICTION_MAP)$/){
      my $cId=processAlignment($reader);
      $cId=~s/\s+.*$//; # remove anything after the space
      $contig{$cId}++;
      $reader->next;
    }
  }

  # well placed contigs are those that are mentioned one time as a restriction map
  # and optionally one time as they are mapped to the in situ restriction map
  # meaning, 3 or more mentions and it is misassembled
  my $placedContigs;
  while(my($contig,$count)=each(%contig)){
    $$placedContigs{$contig}=1 if($count<3);
  }
  return $placedContigs;
}

sub processAlignment{
  my($reader)=@_;
  #my $reader=XML::LibXML::Reader->new( string => $xml);

  if($reader->name=~/RESTRICTION_MAP/){
    while($reader->read){
      # this line could be buggy because it might also be the restriction_map ID
      if($reader->name=~/SEQUENCE_ANNOTATION/){
        my $id=$reader->getAttribute('IDENTIFIER');
        $id=0 if($reader->getAttribute('INSILICO') =~ /false/i);
        $reader->read while($reader->name !~/RESTRICTION_MAP/);
        return $id;
      }
    }
  }
  elsif($reader->name=~/MAP_ALIGNMENT/){
    my ($id1,$id2)=($reader->getAttribute('MAP1'),$reader->getAttribute('MAP2'));
    for my $id($id1,$id2){
      return $id if($id=~/in silico/);
    }
  }

  die "Error: Could not find the correct contig name. Dump of XML below:\n".$reader->readOuterXml;
  return 0;
}

sub usage{
  "incorporate an OpGen xml to make an assembly better: remove misassembled contigs and perform reference assembly of unused reads back onto the well-assembled contigs
  Usage: perl $0 -r restrictionMap.xml -s sequenceReads.fasta -a assembly.ace -o assembly.fasta
  -r restrictionMap.xml
    Restriction map xml outputted by MapSolver
  -a assembly.ace
    The ace file that shows where each read maps to in the input assembly
  -o assembly.fasta
    The output assembly
  -s sequenceReads.fasta
    The raw reads file in fastq format. If .fasta.qual file exists, it will be used

  Optional
  -t tmp
    The location of the temporary directory
  -k
    Keep the temp files around
  -h
    This help menu
  ";
}
