#!/usr/bin/env perl

# run-assembly: Perform standard assembly protocol operations on 454 pyrosequenced flowgram file(s)
# Author: Andrey Kislyuk (kislyuk@gatech.edu)

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
	$$settings{min_out_contig_length} = 100; #TODO put this into config file

	die("Usage: $0 input.fastq [, input2.fastq, ...] [-R references.mfa] [-C workdir]\n  Input files can also be fasta files.  *.fasta.qual files are considered as the relevant quality files") if @ARGV < 1;

	my @cmd_options = ('ChangeDir=s', 'Reference=s@', 'keep', 'tempdir=s', 'output=s');
	GetOptions($settings, @cmd_options) or die;

	$$settings{outfile} = $$settings{output};
	$$settings{outfile} ||= "$0.out.fasta";
	$$settings{outfile} = File::Spec->rel2abs($$settings{outfile});
	open(FH, '>', $$settings{outfile}) or die("Error writing to output file $$settings{outfile}");
	close FH;
	logmsg "Output file is " . $$settings{outfile} . "\n";

	my @ref_files = @{$$settings{Reference}} if defined $$settings{Reference};
	logmsg "No reference files supplied. Reverting to assembly mode" unless @ref_files;

	my @input_files = @ARGV;
  die "Error: No input files given" if(@ARGV<1);

	$$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($$settings{keep}));
	logmsg "Temporary directory is $$settings{tempdir}";

	foreach my $file (@input_files, @ref_files) {
		$file = File::Spec->rel2abs($file);
		die("Input or reference file $file not found") unless -f $file;
	}

  my $fastqfiles=baseCall(\@input_files,$settings); #array of [0=>fna,1=>qual file] (array of arrays)

	my $final_seqs;

	if (@ref_files) {
    # TODO use Columbus
    die "Using references for illumina assembly is not working right now";
	} else {
    my($velvet_basename,$velvet_assembly,$combined_filename);

    # TODO multithread assemblies

    $velvet_basename = runVelvetAssembly($fastqfiles,$settings);
    $velvet_assembly="$velvet_basename/contigs.fa";

    $combined_filename=$velvet_assembly;

		# TODO: reprocess repeat or long singleton reads
    # repeating the process might be moot if we use reconciliator -LK
		$final_seqs = AKUtils::readMfa("$combined_filename");

	}

	foreach my $seq (keys %$final_seqs) {
		delete $$final_seqs{$seq} if length($$final_seqs{$seq}) < $$settings{min_out_contig_length};
	}
  logmsg "";
	AKUtils::printSeqsToFile($final_seqs, $$settings{outfile}, {order_seqs_by_name => 1});

	logmsg "Output is in $$settings{outfile}";

	if ($$settings{keep}) {
		my ($out_filename, $out_dirname) = fileparse($$settings{outfile});
		logmsg "Saving assembly working directory $$settings{tempdir} to $out_dirname";
		move($$settings{tempdir}, $out_dirname) or die "Error moving output directory $$settings{tempdir} to $out_dirname: $!";
	}

	return 0;
}

# determine file types
sub ext($){ # return the extension of a filename
  my($file)=@_;
  my $ext=$file;
  $ext=~s/^.+\.//;
  return lc($ext);
}
sub is_sff($){
  my($file)=@_;
  return 1 if(ext($file) eq 'sff');
  return 0;
}
sub is_fasta($){
  my($file)=@_;
  return 1 if(ext($file)=~/^(sff|fna|fasta|ffn)$/i);
  return 0;
}
sub is_fastq($){
  my($file)=@_;
  return 1 if(ext($file)=~/^fastq|fq$/i);
}
# wrapper for all base calling for input files
sub baseCall($$){
  my($input_files,$settings)=@_;
  my $fastqfiles;
  foreach my $file (@$input_files){
    if(is_fastq($file)){
      push(@$fastqfiles,$file);
    }
    elsif(is_fasta($file)){
      die "I have not implemented fasta files yet (offending file was $file)";
    }
    else{
      die "The format of $file is incompatible";
    }
  }
  return $fastqfiles;
}

# creates qual and sequence fasta files from FASTQ illumina file (basecalling)
sub fastq2fastaqual($$) {
  die "This subroutine is deprecated in favor of using FASTQ files as input";
	my ($fastq_files, $settings) = @_;
	my @fastaqualfiles;
	foreach my $input_file (@$fastq_files) {
    logmsg "Converting $input_file to fasta/qual files via BioPerl";
		my ($fastq_file, $fastq_dir) = fileparse($input_file);
    my $seqin=Bio::SeqIO->new(-file=>"$fastq_dir/$fastq_file",-format=>"fastq-illumina");
    my $seqout=Bio::SeqIO->new(-file=>">$$settings{tempdir}/$fastq_file.fasta",-format=>"fasta");
    my $qualout=Bio::SeqIO->new(-file=>">$$settings{tempdir}/$fastq_file.qual",-format=>"qual");
    while(my $seq=$seqin->next_seq){
      $seqout->write_seq($seq);
      $qualout->write_seq($seq);
    }
		push(@fastaqualfiles, ["$$settings{tempdir}/$fastq_file.fasta", "$$settings{tempdir}/$fastq_file.qual"]);
	}
	return \@fastaqualfiles;
}

sub runVelvetAssembly($$){
  my($fastqfiles,$settings)=@_;
  my $run_name = "$$settings{tempdir}/velvet";
  system("mkdir -p $run_name") if(!-d $run_name);
  logmsg "Executing Velvet assembly $run_name";
  my $velvetPrefix=File::Spec->abs2rel("$run_name/auto");
  my $command="VelvetOptimiser.pl -a -v -p $velvetPrefix ";
  #$command.="-k 'n50*Lcon' -c 'n50*Lcon' "; # change the metrics for how velvet rates its assembly
  $command.="-f '";
  foreach my $fqFiles (@$fastqfiles){
    my $fileFormat="fastq"; # per the Velvet manual
    #TODO detect the chemistry of each run and treat them differently (454, Illumina, etc)
    #see if the reads are mate pairs (454)
    my $isMatePairs=0;
    # TODO find out if there are mate pairs (observed-insert-lengths program via Velvet?)
    #TODO remove reads with an overall bad quality (about <20)
    #TODO examine reads to see if the file overall belongs in "long" "short" etc: might just filter based on chemistry
    my $readLength="short"; # per velvet docs
    $command.="-$readLength -$fileFormat $fqFiles ";
  }
  $command.="' 2>&1 ";
  logmsg "VELVET COMMAND\n  $command";
  system($command); die if $?;

  # cleanup
  my @velvetTempDir=glob("$velvetPrefix*"); # find the output directory
  logmsg "Velvet directory(ies) found: ".join(", ",@velvetTempDir) ." . Assuming $velvetTempDir[0] is the best Velvet assembly.";
  system("cp -r $velvetTempDir[0]/* $run_name/"); # copy the output directory contents to the actual directory

  #TODO create dummy qual file (phred qual for each base is probably about 60-65). Or, if Velvet outputs it in the future, add in the qual parameter.

  system("amos2ace $run_name/velvet_asm.afg"); logmsg "Warning: amos2ace failed" if $?; # make an ace file too

  return $run_name;
}

# TODO use run_assembly_metrics.pl instead to streamline
sub count_contigs{
	my $file=shift;
	open(FH,"<$file")or die "Could not find $file because $!";
	my @lines=<FH>;
	close(FH);
	my @contigs=grep /^>/,@lines;
	return scalar @contigs;
}
