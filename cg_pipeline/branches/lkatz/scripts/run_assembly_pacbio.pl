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

	my @cmd_options = qw(keep tempdir=s outfile=s illuminaReads=s@ sffReads=s@ ccsReads=s@);
	GetOptions($settings, @cmd_options) or die;

	# STEP 1 receive input files
  $$settings{outfile} ||= "$0.out.fasta";
	$$settings{outfile} = File::Spec->rel2abs($$settings{outfile});
	open(FH, '>', $$settings{outfile}) or die("Error writing to output file $$settings{outfile}");
	close FH;
	logmsg "Output file is " . $$settings{outfile} . "\n";

	my @input_files = @ARGV;
  die "Error: No input files given" if(@ARGV<1);

  # set up the temp directory
	$$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($$settings{keep}));
  mkdir($$settings{tempdir}) if(!-d $$settings{tempdir});
  $$settings{tempdir}=File::Spec->rel2abs($$settings{tempdir});
	logmsg "Temporary directory is $$settings{tempdir}";

  # make all input file abs paths and check that they exist
	foreach my $file (@input_files,@{$$settings{illuminaReads}}, @{$$settings{sffReads}}, @{$$settings{ccsReads}} ) {
		$file = File::Spec->rel2abs($file);
		die("Input file $file not found") unless -f $file;
	}
  
  # STEP 2 modify input files and create FRG files
  my $inputDir=highQualityReadsToFrg(\@input_files,$settings);

  # STEP 3 create FRG file for long reads
  padLongReads($inputDir,$settings);

  # STEP 4 assembly
  my $caAssembly=caAssembly($inputDir,$settings);

  # STEP 5 AHA?

  my $final_seqs = AKUtils::readMfa($caAssembly);

	foreach my $seq (keys %$final_seqs) {
		delete $$final_seqs{$seq} if length($$final_seqs{$seq}) < $$settings{assembly_min_contig_length};
	}
	AKUtils::printSeqsToFile($final_seqs, $$settings{outfile}, {order_seqs_by_name => 1});

	logmsg "Output is in $$settings{outfile}";

	return 0;
}

sub highQualityReadsToFrg{
  my($input_files,$settings)=@_;
  my $newInputDir="$$settings{tempdir}/input";
  mkdir($newInputDir);
  # STEP 2a standardize pacbio fastq files
  logmsg "Reading and cleaning long reads";
  my $longreadsFile="$newInputDir/longreads.fastq";
  if(!-e $longreadsFile){
    my $longreadsInputDir="$newInputDir/longreads";
    mkdir("$longreadsInputDir");
    for(my $i=0;$i<@$input_files;$i++){
      my $file=$$input_files[$i];
      my $filename=basename($file);
      my $standardizedFile="$longreadsInputDir/$filename";
      command("run_assembly_convertMultiFastqToStandard.pl -i '$file' -o '$standardizedFile' -m 50 2>&1") if(!-e $standardizedFile);
      $$input_files[$i]=$standardizedFile;
    }

    # instead of just concatenating all reads, see if they can be filtered a little first
    my $command="run_assembly_trimClean.pl --min_avg_quality 6 --notrim --min_length 50 ";
    $command.="-i '$_' " for(@$input_files);
    $command.="-o $longreadsFile 2>&1 ";
    command($command);
    # Remove these temporary intermediate files.  They are huge and not necessary, even in a temp directory
    command("rm -rf $longreadsInputDir");
  }

  # STEP 2b remove CCS reads from long reads (ie remove duplication)
  logmsg "Reading and cleaning CCS pacbio reads";
  my $ccsreadsFile="$newInputDir/ccsreads.fastq";
  if(!-e $ccsreadsFile){
    my $command="run_assembly_trimClean.pl --min_avg_quality 20 --notrim --min_length 50 ";
    $command.="-i '$_' " for(@{$$settings{ccsReads}});
    $command.="-o $ccsreadsFile 2>&1";
    command($command);
  }


  # STEP 2c create FRG files for each input file except pacbio long reads.
  # All FRG files found in the input directory will be assumed to be part
  # of the assembly.  Do not contaminate the input directory!
  
  # CCS reads
  for my $fastq (@{$$settings{ccsReads}}){
    ccsToFrg($fastq,$newInputDir,$settings);
  }
  # illumina
  for my $fastq (@{$$settings{illuminaReads}}){
    my $frg=illuminaToFrg($fastq,$newInputDir,$settings);
  }
  # 454
    # sffToCA
  warn("WARNING: sffToCA has not been implemented");

  die "Done with step 2c. I will not continue";

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

##########
## Utility
##########

sub ccsToFrg{
  my($fastq,$newInputDir,$settings)=@_;
  my $libraryname=basename($fastq,qw(.fastq .fq));
  my $frgFile="$newInputDir/$libraryname.frg";
  if(-e $frgFile){
    warn "WARNING: frg file $frgFile already exists. I will not overwrite it";
    return $frgFile;
  }

  # fastq to fasta/qual
  logmsg "Creating a fasta/qual for $fastq";
  my ($fasta,$qual)=("$newInputDir/$libraryname.fna","$newInputDir/$libraryname.qual");
  open(IN,$fastq) or die "Could not open $fastq: $!";
  open(FASTA,">$fasta") or die "Could not open $fasta for writing: $!";
  open(QUAL,">$qual") or die "Could not open $qual for writing: $!";
  while(my $id=<IN>){
    $id=~s/^@//;
    my $sequence=<IN>;
    <IN>; # discard the + line
    my $qual=<IN>;
    chomp($qual);
    my @qual=split(//,$qual);
    @qual=map(ord($_)-33,@qual);
    $qual=join(" ",@qual);
    
    print FASTA ">$id$sequence";
    print QUAL  ">$id$qual\n";
  }
  close QUAL;
  close FASTA;
  close IN;

  # fasta/qual to frg
  my $command="convert-fasta-to-v2.pl -pacbio -s '$fasta' -q '$qual' -l $libraryname > '$frgFile'";
  command($command);

  return $frgFile;
}
sub illuminaToFrg{
  my($fastq,$newInputDir,$settings)=@_;
  my $libraryname=basename($fastq,qw(.fastq .fq));
  my $frgFile="$newInputDir/$libraryname.illumina.frg";
  if(-e $frgFile){
    warn "WARNING: frg file $frgFile already exists. I will not overwrite it";
    return $frgFile;
  }

  my $is_PE=is_illuminaPE($fastq,$settings);
  my $command="fastqToCA -insertsize 300 20 -libraryname $libraryname -technology illumina -innie ";
  if($is_PE){
    $command.="-mates '$fastq' ";
  } else {
    $command.="-reads '$fastq' ";
  }
  $command.="1> $frgFile";
  command($command);
  return $frgFile;
}

# return whether or not the input file is a shuffled PE fastq
sub is_illuminaPE{
  my($fastq,$settings)=@_;
  # 20 reads is probably enough to make sure that it's shuffled
  my $numEntriesToCheck=$$settings{pefastq_numEntriesToCheck}||20;
  my $numEntries=0;
  open(IN,$fastq) or die "Could not open $fastq for reading: $!";
  while(my $read1Id=<IN>){
    my $discard;
    $discard=<IN> for(1..3);
    my $read2Id=<IN>;
    $discard=<IN> for(1..3);

    if($read1Id!~/\/1$/ || $read2Id!~/\/2$/){
      return 0;
    }

    $numEntries++;
    last if($numEntries>=$numEntriesToCheck);
  }
  close IN;

  return 1;
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
