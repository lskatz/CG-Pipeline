#!/usr/bin/env perl

# run-assembly-filterContigs: remove contigs that fit certain criteria
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
use File::Temp ('tempdir');
use File::Path;
use File::Spec;
use File::Basename;
use Data::Dumper;

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {print STDERR "$0: @_\n";}

exit(main());

sub main() {
	$settings = AKUtils::loadConfig($settings);

	my @cmd_options = qw(verbose assembly_min_contig_length|length=i help numcpus=i rename);
	GetOptions($settings, @cmd_options) or die;
	$$settings{numcpus}||=1;
  $$settings{assembly_min_contig_length}||=500;
	die usage($settings) if @ARGV < 1;
	die usage($settings) if ($$settings{help});

  my %seq;
  for(@ARGV){
    my $seq=AKUtils::readMfa($_,$settings);
    $seq=filterSeqs($seq,$settings);
    %seq=(%seq,%$seq);
  }

  while(my($id,$sequence)=each(%seq)){
    $sequence=~s/(.{80})/$1\n/g;
    $sequence=~s/^\s+|\s+//g;     # trim the sequence's whitespace
    print ">$id\n$sequence\n";    # but put a newline after the sequence
  }

  return 0;
}

sub filterSeqs{
  my($seq,$settings)=@_;
  my %seq;
  my $i=0;
  while(my($id,$sequence)=each(%$seq)){
    if(length($sequence) < $$settings{assembly_min_contig_length}){
      logmsg "Removing short contig $id" if($$settings{verbose});
      next;
    }
    
    # Get rid of uncomplex contigs:
    #   Contigs of a K-mer repeat
    my $uncomplexkmer="";
    for(my $k=1;$k<20;$k++){
      my $kmer=substr($sequence,0,$k);
      if($sequence=~/^($kmer)+$/i){
        $uncomplexkmer=$kmer;
        last;
      }
    }
    if($uncomplexkmer){
      logmsg "Removing uncomplex contig $id with kmer $uncomplexkmer";
      next;
    }

    if($$settings{rename}){
      my $newid = sprintf("contig%06d",++$i);
      $id=$newid;
      logmsg "Renaming $id to $newid" if($$settings{verbose});
    }

    $seq{$id}=$sequence;
  }
  return \%seq;
}

sub usage{
  my ($settings)=@_;
  "Filters sequences from a fasta file that do not meet certain criteria and prints them to stdout
  Usage: $0 file.fasta [file2.fasta...] > out.fasta
  --length  $$settings{assembly_min_contig_length}   Minimum size of a contig
  --rename        rename contigs
  --verbose       
  "
}

