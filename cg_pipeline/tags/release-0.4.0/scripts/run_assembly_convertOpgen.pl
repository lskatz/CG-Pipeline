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
use Bio::Seq::Quality;
use Math::Round qw(round nearest);

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
  $$settings{insertIncrement}||=10; # in bp, the distance from one read to the next
  # the linker sequence was obtained from
  # http://seqanswers.com/forums/showthread.php?t=12940
  $$settings{insertSequence}||="GTTGGAACCGAAAGGGTTTGAATTCAAACCCTTTCGGTTCCAAC";

  my @restrictionPos=readMap($restrictionMapDocument,$settings);
  
  # make a raw reads file
  my $reads=writeFastqReads($restrictionSequence,\@restrictionPos,$settings);
  #my $reads=write454Reads($restrictionSequence,\@restrictionPos,$settings);

  move($reads,$outfile);
  logmsg "reads are in $outfile";

  return 0;
}

sub writeSangerReads{
  #http://contig.wordpress.com/2011/01/21/newbler-input-ii-sequencing-reads-from-other-platforms/
}

sub write454Reads{
  my($restrictionSequence,$rPos,$settings)=@_;
  my $dir=$$settings{tempdir};
  my (%outseq,%outqual);
  
  my $sequence=join("",$restrictionSequence,$$settings{insertSequence},$restrictionSequence);
  my $restrictionQual="$$settings{qualScore} " x length($sequence);
  chomp($restrictionQual);
  
  # Make a forward and reverse read.
  # Put in a special linker insert (need to figure what that is)
  for(my $c=0;$c<@$rPos;$c++){
    my $d=$c+1;
    my $rPos=$$rPos[$c];
    for(my $p=0;$p<@$rPos;$p++){
      my $insLength=nearest($$settings{insertIncrement},$$rPos[$p]);
      if(!$outseq{$insLength}){
        $outseq{$insLength}=Bio::SeqIO->new(-format=>"fasta",-file=>">$dir/$insLength.fasta");
        $outqual{$insLength}=Bio::SeqIO->new(-format=>"qual",-file=>">$dir/$insLength.qual");
      }
      my $id=join("_","read",$d,$p);
      my $seq=Bio::Seq::Quality->new(-id=>$id,-seq=>$sequence,-qual=>$restrictionQual);
      #my $qual=Bio::Seq::Quality->new(-id=>$id,-qual=>$restrictionQual);
      $outseq{$insLength}->write_seq($seq);
      $outqual{$insLength}->write_seq($seq);
    }
  }

  return $dir;
}

sub writeFastqReads{
  my($restrictionSequence,$rPos,$settings)=@_;

  my $dir=$$settings{tempdir};
  my %outseq;
  my $restrictionQual="$$settings{qualScore} " x length($restrictionSequence);
  chomp($restrictionQual);
  my $machineIdentifier=$$settings{appname};

  # Make a forward read for each position, matching with the next position up.
  # Make a reverse read for each position, matching with the previous position.
  # Transitively, there will be reads going back to this site too.
  for(my $c=0;$c<@$rPos;$c++){
    my $d=$c+1; # for sequence channel ID
    my $rPos=$$rPos[$c];
    for(my $p=0;$p<@$rPos;$p++){
      my $insLength=nearest($$settings{insertIncrement},$$rPos[$p]);
      if(!$outseq{$insLength}){
        $outseq{$insLength}=Bio::SeqIO->new(-format=>"fastq-illumina",-file=>">$dir/$insLength.fastq");
      }
      my $fid=join("",$machineIdentifier,$p,"#$d",'/1');
      my $rid=join("",$machineIdentifier,$p,"#$d",'/2');
      my $fSeq=Bio::Seq::Quality->new(-id=>$fid,-seq=>$restrictionSequence,-qual=>$restrictionQual);
      my $rSeq=Bio::Seq::Quality->new(-id=>$rid,-seq=>$restrictionSequence,-qual=>$restrictionQual);
      $outseq{$insLength}->write_seq($fSeq,$rSeq);
    }
  }

  logmsg "Created ".scalar(keys(%outseq))." paired end files";
    
  return $dir;
}

sub readMap{
  my($map,$settings)=@_;
 
  my @rLength;

  my $reader = new XML::LibXML::Reader(location => $map) or die "cannot read $map\n";
  while($reader->read){
    if($reader->name eq 'RESTRICTION_MAP' && $reader->getAttribute('INSILICO') eq 'false'){
      my $chrRMap=[];
      while($reader->read){
        if($reader->name eq 'RESTRICTION_MAP'){ # break out at the closing tag
          push(@rLength,$chrRMap);
          last;
        }
        next if($reader->name ne 'F'); # F for fragment
        push(@$chrRMap,($reader->getAttribute('S')+0));
      }
    }
  }

  my @rLocation;
  my $siteCounter=0;
  my $c=0; # chromosome counter
  for my $chrRMap(@rLength){
    my $chrLocation=[0];
    for(my $r=1;$r<@$chrRMap;$r++){
      $$chrLocation[$r]=$rLength[$c][$r]+$$chrLocation[$r-1];
      $siteCounter++;
    }
    $rLocation[$c++]=$chrLocation;
  }
  logmsg "Found $siteCounter sites";

  logmsg "Returning lengths";return @rLength;
  return @rLocation;
}

# debugs an XML node
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

