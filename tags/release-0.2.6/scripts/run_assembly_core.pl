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

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH = $FSFind::LOG || *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

exit(main());

sub main() {
	$settings = AKUtils::loadConfig($settings);

	die("Usage: $0 input.sff [, input2.sff, ...] [-R references.mfa] [-C workdir]") if @ARGV < 1;

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

	$$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($$settings{keep}));
	logmsg "Temporary directory is $$settings{tempdir}";

	foreach my $file (@input_files, @ref_files) {
		$file = File::Spec->rel2abs($file);
		die("Input or reference file $file not found") unless -f $file;
	}

	my $final_seqs;
	if (@ref_files) {
		my $newbler_basename = runNewblerMapping(\@input_files, \@ref_files, $settings);

		my $fastaqualfiles = sff2fastaqual(\@input_files, $settings);

		my $amos_basename = runAMOSMapping($fastaqualfiles, \@ref_files, $settings);

		my $combined_filename = combineNewblerAMOS($newbler_basename, $amos_basename, $settings);

		$final_seqs = AKUtils::readMfa($combined_filename);
	} else {
		my $newbler_basename = runNewblerAssembly(\@input_files, $settings);

		# TODO: reprocess repeat or long singleton reads
		$final_seqs = AKUtils::readMfa("$newbler_basename/assembly/454AllContigs.fna");
	}
	
	$$settings{min_out_contig_length} = 500;

	foreach my $seq (keys %$final_seqs) {
		delete $$final_seqs{$seq} if length($$final_seqs{$seq}) < $$settings{min_out_contig_length};
	}
	AKUtils::printSeqsToFile($final_seqs, $$settings{outfile}, {order_seqs_by_name => 1});

	#copy($combined_filename, $$settings{outfile}) or die("Error writing to output file $$settings{outfile}");
	logmsg "Output is in $$settings{outfile}";

	if ($$settings{keep}) {
		my ($out_filename, $out_dirname) = fileparse($$settings{outfile});
		logmsg "Saving assembly working directory $$settings{tempdir} to $out_dirname";
		move($$settings{tempdir}, $out_dirname) or die "Error moving output directory $$settings{tempdir} to $out_dirname: $!";
	}

	return 0;
}

# creates qual and sequence fasta files for an SFF file
sub sff2fastaqual($$) {
	my ($sff_files, $settings) = @_;
	my @fastaqualfiles;
	foreach my $input_file (@$sff_files) {
		my ($sff_file, $sff_dir) = fileparse($input_file);
		my $invoke_string = "sffinfo -seq '$sff_dir/$sff_file' > '$$settings{tempdir}/$sff_file.fasta'";
		logmsg "Running $invoke_string";
		system($invoke_string); die if $?;
		$invoke_string = "sffinfo -qual '$sff_dir/$sff_file' > '$$settings{tempdir}/$sff_file.qual'";
		logmsg "Running $invoke_string";
		system($invoke_string); die if $?;
		push(@fastaqualfiles, ["$$settings{tempdir}/$sff_file.fasta", "$$settings{tempdir}/$sff_file.qual"]);
	}
	return \@fastaqualfiles;
}

sub runNewblerMapping($$$) {
	my ($input_files, $ref_files, $settings) = @_;
	my $run_name = "$$settings{tempdir}/P__runMapping";
	logmsg "Executing Newbler mapping project $run_name";

	system("newMapping '$run_name'"); die if $?;
	foreach my $ref_file (@$ref_files) {
		system("setRef '$run_name' '$ref_file'"); die if $?;
	}
	foreach my $sff_file (@$input_files) {
		system("addRun '$run_name' '$sff_file'"); die if $?;
	}
#sed -i -e 's|\(<overlapMinMatchLength>\)40\(</overlapMinMatchLength>\)|\125\2|' \
#               -e 's|\(<overlapMinMatchIdentity>\)90\(</overlapMinMatchIdentity>\)|\175\2|' \
#               $runName/mapping/454MappingProject.xml
	system("runProject '$run_name'"); die if $?;
	return $run_name;
}

sub runNewblerAssembly($$) {
	my ($input_files, $settings) = @_;
	my $run_name = "$$settings{tempdir}/P__runAssembly";
	logmsg "Executing Newbler assembly project $run_name";

	system("newAssembly '$run_name'"); die if $?;
	foreach my $sff_file (@$input_files) {
		system("addRun '$run_name' '$sff_file'"); die if $?;
	}
#sed -i -e 's|\(<overlapMinMatchLength>\)40\(</overlapMinMatchLength>\)|\125\2|' \
#               -e 's|\(<overlapMinMatchIdentity>\)90\(</overlapMinMatchIdentity>\)|\175\2|' \
#               $runName/mapping/454MappingProject.xml
	system("runProject '$run_name'"); die if $?;
	return $run_name;
}

sub runAMOSMapping($$$) {
	my ($input_files, $ref_files, $settings) = @_;
	my $run_name = "$$settings{tempdir}/amos_mapping";
	logmsg "Executing AMOS mapping project $run_name";

	my @afg_files;
	foreach my $input_file_pair (@$input_files) {
		my ($input_fasta_file, $input_qual_file) = @$input_file_pair;
		system("toAmos -s '$input_fasta_file' -q '$input_qual_file' -o '$input_fasta_file.afg'"); die if $?;
		push(@afg_files, "$input_fasta_file.afg");
	}

	my $invoke_string = "AMOScmp";
	$invoke_string .= " -D TGT='$_'" for @afg_files;
	$invoke_string .= " -D REF='$_'" for @$ref_files;
	$invoke_string .= " $run_name";
	logmsg "Running $invoke_string";
	system($invoke_string); die if $?;

	return $run_name;
}

sub combineNewblerAMOS($$$) {
	my ($newbler_basename, $amos_basename, $settings) = @_;
	logmsg "Running Newbler-AMOS combining on $newbler_basename, $amos_basename";

	my (%unmapped_reads, %repeat_reads);
	open(IN, '<', "$newbler_basename/mapping/454ReadStatus.txt") or die;
	while (<IN>) {
		chomp;
		my ($read_id, $status) = split /\s+/;
		$unmapped_reads{$read_id} = 1 if $status eq 'Unmapped';
		$repeat_reads{$read_id} = 1 if $status eq 'Repeat';
	}
	close IN;

	open(OUT, '>', "$$settings{tempdir}/newbler_unmapped_acc") or die;
	print OUT "$_\n" for keys %unmapped_reads;
	close OUT;
	open(OUT, '>', "$$settings{tempdir}/newbler_repeat_acc") or die;
	print OUT "$_\n" for keys %repeat_reads;
	close OUT;

	system("sfffile -i '$$settings{tempdir}/newbler_unmapped_acc' -o '$$settings{tempdir}/newbler_unmapped.sff' '$newbler_basename'/sff/*.sff");
	die if $?;
	my $unmapped_fastaqual = sff2fastaqual(["$$settings{tempdir}/newbler_unmapped.sff"], $settings);
	system("sfffile -i '$$settings{tempdir}/newbler_repeat_acc' -o '$$settings{tempdir}/newbler_repeat.sff' '$newbler_basename'/sff/*.sff");
	die if $?;
	my $repeat_fastaqual = sff2fastaqual(["$$settings{tempdir}/newbler_repeat.sff"], $settings);

	my $combined_fasta_file = "$$settings{tempdir}/combined_in.fasta";
	system("cat '$newbler_basename/mapping/454AllContigs.fna' '$amos_basename.fasta' '$$unmapped_fastaqual[0]->[0]' > $combined_fasta_file");
#	system("cat '$newbler_basename/mapping/454AllContigs.fna' '$amos_basename.fasta' > $combined_fasta_file");
	die if $?;
	system("toAmos -s '$combined_fasta_file' -o '$combined_fasta_file.afg'");
	die if $?;
	system("minimus -D TGT='$combined_fasta_file.afg' '$$settings{tempdir}/minimus.combined'");
	die if $?;

	return "$$settings{tempdir}/minimus.combined.fasta";
}
