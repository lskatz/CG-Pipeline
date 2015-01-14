#!/usr/bin/env perl
# Wrapper script for all cleaning reads.

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
use AKUtils;

use Getopt::Long;
use File::Temp ('tempdir');
use File::Path;
use File::Spec;
use File::Copy;
use File::Basename;
use List::Util qw(min max sum shuffle);
use Data::Dumper;

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH=*STDERR; print $FH "$0: ".(caller(1))[3].": @_\n";}

exit(main());

sub main() {
  GetOptions($settings, qw(help tempdir=s trimCleanXopts=s strimUpAndDownXopts=s duplicatesXopts=s));
  die usage() if($$settings{help});
  $$settings{tempdir}||=AKUtils::mktempdir();
  $$settings{trimCleanXopts}||="";
  $$settings{strimUpAndDownXopts}||="";
  $$settings{duplicatesXopts}||="";
  
  mkdir $$settings{tempdir} if(!-d $$settings{tempdir});
  my @read=@ARGV;
  die usage() if(!@read);

  logmsg "Trimming all files into a single comprehensive file";
  my $trimmedReads=trimUpAndDown(\@read,$settings);
  logmsg "Removing duplicate reads";
  my $noDupes=removeDuplicateReads($trimmedReads,$settings);
  logmsg "Cleaning reads";
  my $cleaned=trimClean($noDupes,$settings);

  system("cat '$cleaned'"); die if $?;

  return 0;
}

sub trimUpAndDown{
  my($read,$settings)=@_;
  my $trimmedReads="$$settings{tempdir}/reads.trimmed.fastq.gz";
  system("run_assembly_trimLowQualEnds.pl $$settings{strimUpAndDownXopts} ".join(" ",@$read)." | gzip -c > $trimmedReads");
  die if $?;
  return $trimmedReads;
}

sub removeDuplicateReads{
  my($reads,$settings)=@_;
  my $noDupes="$$settings{tempdir}/reads.noDupes.fastq.gz";
  system("run_assembly_removeDuplicateReads.pl $$settings{duplicatesXopts} $reads > $noDupes");
  die if $?;
  return $reads;
}

sub trimClean{
  my($reads,$settings)=@_;
  my $out="$$settings{tempdir}/reads.trimmedCleaned.fastq";
  system("run_assembly_trimClean.pl $$settings{trimCleanXopts} -i $reads -o $out 1>&2 ");
  die if $?;
  return $out;
}

sub usage{
  local $0=basename $0;
  "This is a wrapper around all CGP cleaning scripts.
  Usage: $0 reads.fastq [reads2.fastq...] > reads.cleaned.fastq
    -h for this help menu

    The following options can send options directly to the sub-scripts. 
    These parameters are not sanity-checked.
    -t ' ' run_assembly_trimClean.pl options
    -s ' ' run_assembly_trimLowQualEnds.pl options
    -d ' ' run_assembly_removeDuplicateReads.pl options
  "
}
