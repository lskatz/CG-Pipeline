#!/usr/bin/env perl

# run-assembly-pacbio: Perform standard assembly protocol operations on pacbio files
# Author: Lee Katz <lkatz@cdc.gov>
# Author: Chandni Desai <CDesai@cdc.gov>
# Author: Satish Ganakammal <Satishkumar_ranganathanganakammal@sra.com>

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
use Bio::Perl;
use Data::Dumper;

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH = $FSFind::LOG || *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

exit(main());

sub main() {
	$settings = AKUtils::loadConfig($settings);
	$$settings{assembly_min_contig_length} ||= 500; #TODO put this into config file

	die(usage()) if @ARGV < 1;

	my @cmd_options = qw(keep tempdir=s outfile=s illuminaReads=s@ sffReads=@ ccsReads=s@);
	GetOptions($settings, @cmd_options) or die;

	# STEP 1 receive input files
  $$settings{outfile} ||= "$0.out.fasta";
	$$settings{outfile} = File::Spec->rel2abs($$settings{outfile});
	open(FH, '>', $$settings{outfile}) or die("Error writing to output file $$settings{outfile}");
	close FH;
	logmsg "Output file is " . $$settings{outfile} . "\n";

	my @input_files = @ARGV;
  die "Error: No input files given" if(@ARGV<1);

	$$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($$settings{keep}));
  mkdir($$settings{tempdir}) if(!-d $$settings{tempdir});
  $$settings{tempdir}=File::Spec->rel2abs($$settings{tempdir});
	logmsg "Temporary directory is $$settings{tempdir}";

	foreach my $file (@input_files) {
		$file = File::Spec->rel2abs($file);
		die("Input or reference file $file not found") unless -f $file;
	}
  
  # STEP 2 modify input files
  my $inputDir=modifyInputFiles(\@input_files,$settings);

  # STEP 3 create FRG file for long reads
  padLongReads($inputDir,$settings);

  # STEP 4 assembly
  caAssembly($inputDir,$settings);

  # STEP 5 AHA?

  $final_seqs = AKUtils::readMfa("$combined_filename");

	foreach my $seq (keys %$final_seqs) {
		delete $$final_seqs{$seq} if length($$final_seqs{$seq}) < $$settings{assembly_min_contig_length};
	}
	AKUtils::printSeqsToFile($final_seqs, $$settings{outfile}, {order_seqs_by_name => 1});

	logmsg "Output is in $$settings{outfile}";

	return 0;
}

sub modifyInputFiles{
  my($input_files,$settings)=@_;
  my $newInputDir="$$settings{tempdir}/input";
  mkdir($newInputDir);
  # STEP 2a standardize pacbio fastq files
  my $longreadsInputDir="$newInputDir/longreads";
  mkdir("$longreadsInputDir");
  for(my $i=0;$i<@$input_files;$i++){
    my $file=$$input_files[$i];
    my $filename=basename($file);
    my $standardizedFile="$longreadsInputDir/$filename";
    command("run_assembly_convertMultiFastqToStandard.pl -i '$file' -o '$standardizedFile' 2>&1") if(!-e $standardizedFile);
    $$input_files[$i]=$standardizedFile;
  }
  # STEP 2b remove CCS reads from long reads (ie remove duplication)
  # STEP 2c create FRG files for each input file except pacbio long reads
  # TODO use logic to find if a file is a shuffled paired end file

  return $newInputDir;
}

sub padLongReads{
  my($inputDir,$settings)=@_;
  # STEP 3a create or locate appropriate spec files
  # STEP 3b pad the long reads with pacBioToCA
}
sub caAssembly{
  my($inputDir,$settings);
  # STEP 4a create or locate appropriate spec files
  # STEP 4b runCA
}

# run a command
sub command{
  my ($command)=@_;
  logmsg "RUNNING COMMAND\n  $command";
  system($command);
  die "ERROR running command $command\n  With error code $?. Reason: $!" if($?);
  return 1;
}

sub usage{
	"Usage: $0 input.fastq [, input2.fastq, ...] [-o outfile.fasta] -c pacbio.ccs.fastq -i illumina.fastq -s 454.sff
  Input files should be the filtered subreads from pacbio
  -c, -i, -s
    You can use multiple files. 
    Designate each new file with a new -c, -i, or -s flag.
    Multiple files from any platform are allowed.
    Paired end illumina data can be inputted using shuffled sequences in a singled file.
  -t tempdir
    Designate a different temporary directory
  -k
    If using default tempdir, use -k to keep temporary files
  "
}
