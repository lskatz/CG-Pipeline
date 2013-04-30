#!/usr/bin/env perl

# run_assembly_contigsToScaffolds.pl: combine contigs into scaffolds in a user-specified manner
# Author: Lee Katz (gzu2@cdc.gov)

package PipelineRunner;
my ($VERSION) = ('$Id: $' =~ /,v\s+(\d+\S+)/o);

my $settings = {
    appname => 'cgpipeline',
};

use strict;
no strict "refs";
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
use Bio::SeqUtils;

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH = $FSFind::LOG || *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

exit(main());

sub main() {
  $settings = AKUtils::loadConfig($settings);
  die("Incorrect arguments. Usage is as follows.
  Usage: $0 -i input.fasta -o output.fasta -s contig1:distanceInBp:contig5 -s contig00045:1000:contig00002 --defaultSpacing=10
  Place a negative number in front of the contig to invert, e.g., -contig1:10:-contig2, where both contigs are inverted and the distance is 10 base pairs.
  Default spacing will be used when presented with a double colon (::) and is 10 bp when not provided by the user.
  Scaffold ordering is performed in order, as given by the ordering of the -s arguments.
  ") 
  if @ARGV<1;

  my @cmd_options=('output=s','input=s','scaffold=s@','defaultSpacing=s');
  GetOptions($settings, @cmd_options) or die;

  $$settings{defaultSpacing} ||= 10;
  my %seqHash=seqHash($$settings{input});
  
  my %scaffold=makeScaffoldsFile(\%seqHash,$settings);
  
  logmsg "Output file is $$settings{output}";
  
  return 0;
}

sub makeScaffoldsFile{
  my($seqHash,$settings)=@_;
  
  # make all scaffolds
  my $seqout=makeScaffolds($seqHash,$settings);
  
  # TODO print out the rest to seqout
  printRestOfSeqs($seqHash,$seqout);
}

sub printRestOfSeqs{
  my($seqHash,$seqout)=@_;
  foreach my $seq (values(%$seqHash)){
    $seqout->write_seq($seq);
  }
  
  return 1;
}

sub makeScaffolds{
  my($seqHash,$settings)=@_;
  my $seqout=Bio::SeqIO->new(-file=>">$$settings{output}");
  
  my $scaffoldId=0;
  foreach my $cigarLine (@{ $$settings{scaffold} }){
    my @field=split(/:/,$cigarLine);
    my @seqToConcat=();
    my @id=();
    for(my $i=0;$i<@field;$i+=2){
      my $contigId=$field[$i];
      my $strand='+';
      if($contigId=~/([+-])(.+)/){
        $contigId=$2;
        $strand=$1;
      }
      my $seq=$$seqHash{$contigId} or die "Could not locate the contig $contigId in the input file $$settings{input}";
      $seq=$seq->revcom if($strand eq '-' || $strand==0);
      
      # separator
      my $separator=$field[$i+1]||$$settings{defaultSpacing}; # spacing is what's listed, the default settings, or nothing
      #$separator=0 if($i+1==@field);
      $separator="N"x$separator;
      my $separatorSeq=Bio::Seq->new(-id=>"Separator",-seq=>$separator,-sequence=>$separator);
      
      # mark for concatenation
      push(@seqToConcat,$seq);
      push(@seqToConcat,$separatorSeq) unless($i+1==@field);
      delete($$seqHash{$contigId});
      
      push(@id,"$strand$contigId");
    }
    Bio::SeqUtils->cat(@seqToConcat);
    my $seq=$seqToConcat[0];
    $seq->id("Scaffold".++$scaffoldId."_".join("_",@id));
    $seq->desc(" ");
    $seqout->write_seq($seqToConcat[0]);   
  }
  return $seqout;
}

sub seqHash{
  my($in)=@_;
  my $seqin=Bio::SeqIO->new(-file=>$$settings{input});
  my %hash;
  while(my $seq=$seqin->next_seq){
    $hash{$seq->id}=$seq;
  }
  return %hash;
}