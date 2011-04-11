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
use CGPipeline::Reconciliator;
use Data::Dumper;

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH = $FSFind::LOG || *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

exit(main());

sub main() {
	$settings = AKUtils::loadConfig($settings);

	die("Usage: $0 --readsfilename=path/to/reads.fasta --newblerDir someDirectory1 --newblerDir someDirectory2 --velvetDir someDirectory3 --output=outputDirectory [--tempdir temporaryDirectory] [--force]") if @ARGV < 1;

	my @cmd_options = ('ChangeDir=s', , 'keep', 'tempdir=s', 'output=s', 'newblerDir=s@', 'velvetDir=s@', 'amoscmpDir=s@', 'readsfilename=s','force');
	GetOptions($settings, @cmd_options) or die;

	$$settings{outfile} = $$settings{output};
	$$settings{outfile} ||= "$0.out.fasta";
	$$settings{outfile} = File::Spec->rel2abs($$settings{outfile});
	open(FH, '>', $$settings{outfile}) or die("Error writing to output file $$settings{outfile}");
	close FH;
	logmsg "Output file is " . $$settings{outfile};

	$$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($$settings{keep}));
	logmsg "Temporary directory is $$settings{tempdir}";

  # add all assemblies to the reconciliator
  # TODO do not add an assembly if it already has been converted, unless there is the force command
  my $reconciliator=CGPipeline::Reconciliator->new({outputDir=>$$settings{tempdir}});
  for my $amoscmpDir (@{ $$settings{amoscmpDir} }){
    $reconciliator->addAmoscmpAssembly({amoscmpDir=>$amoscmpDir,readsfasta=>$$settings{readsfilename},force=>$$settings{force}});
  }
  for my $velvetDir (@{ $$settings{velvetDir} }){
    $reconciliator->addVelvetAssembly({velvetDir=>$velvetDir,force=>$$settings{force}});
  }
  for my $newblerDir (@{ $$settings{newblerDir} }){
    $reconciliator->addNewblerAssembly({newblerDir=>$newblerDir,force=>$$settings{force}});
  }
  die "Ending right after reconciliator format conversions\n";
  exit;
  print Dumper $reconciliator;

	#AKUtils::printSeqsToFile($final_seqs, $$settings{outfile}, {order_seqs_by_name => 1});

	#logmsg "Output is in $$settings{outfile}";

	#if ($$settings{keep}) {
	#	my ($out_filename, $out_dirname) = fileparse($$settings{outfile});
#		logmsg "Saving assembly working directory $$settings{tempdir} to $out_dirname";
#		move($$settings{tempdir}, $out_dirname) or die "Error moving output directory $$settings{tempdir} to $out_dirname: $!";
#	}

	return 0;
}

