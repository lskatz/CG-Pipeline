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
sub logmsg {my $FH = $FSFind::LOG || *STDERR; print $FH "$0: ".(caller(1))[3].": @_\n";}

exit(main());

sub main() {
	$settings = AKUtils::loadConfig($settings);

	my @cmd_options = qw(length=i help numcpus=i rename);
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
    if($$settings{rename}){
      $id = sprintf("contig%06d",++$i);
    }
    $seq{$id}=$sequence if(length($sequence) > $$settings{assembly_min_contig_length});
  }
  return \%seq;
}

sub usage{
  my ($settings)=@_;
  "Filters sequences from a fasta file that do not meet certain criteria and prints them to stdout
  Usage: $0 file.fasta [file2.fasta...] > out.fasta
  -l $$settings{assembly_min_contig_length} Minimum size of a contig
  --rename  rename contigs
  "
}
