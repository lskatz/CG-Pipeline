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
use Data::Dumper;

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH = $FSFind::LOG || *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

exit(main());

sub main() {
	$settings = AKUtils::loadConfig($settings);

	die("Usage: $0 input.sff [, input2.sff, ...] [-R references.mfa] [-C workdir]\n  Input files can also be fasta files.  *.fasta.qual files are considered as the relevant quality files") if @ARGV < 1;

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

  my $fastaqualfiles=baseCall(\@input_files,$settings); #array of [0=>fna,1=>qual file] (array of arrays)

	my $final_seqs;

	if (@ref_files) {
    # TODO multithread reference assemblies
		my $newbler_basename = runNewblerMapping(\@input_files, \@ref_files, $settings);

		my $amos_basename = runAMOSMapping($fastaqualfiles, \@ref_files, $settings);

		my $combined_filename = combineNewblerAMOS($newbler_basename, $amos_basename, $settings);

		$final_seqs = AKUtils::readMfa($combined_filename);
	} else {
    
    # TODO multithread assemblies
		my $newbler_basename = runNewblerAssembly(\@input_files, $settings);

    my $combined_filename=combineNewblerDeNovo($newbler_basename,$settings);

		# TODO: reprocess repeat or long singleton reads
    # repeating the process might be moot if we use reconciliator -LK
		$final_seqs = AKUtils::readMfa("$combined_filename");
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
# wrapper for all base calling for input files
sub baseCall($$){
  my($input_files,$settings)=@_;
  my $fastaqualfiles;
  foreach my $file (@$input_files){
    my $seqQualPair; # [0=>seq,1=>qual]
    if(is_sff($file)){
      $seqQualPair=sff2fastaqual([$file], $settings);
      $seqQualPair=$$seqQualPair[0];
    }
    elsif(is_fasta($file)){
      my $qualfile="";
      $qualfile="$file.qual" if(-e "$file.qual");
      $seqQualPair=[$file,$qualfile];
    }
    else{
      die "The format of $file is incompatible";
    }
    push(@$fastaqualfiles,$seqQualPair);
  }
  return $fastaqualfiles;
}

# given a file containing one read accession per line and a pattern for all SFFs, extracts the reads and creates a new SFF
# the pattern for all SFFs is to be used in glob(), e.g. "sff/*.sff"
sub sffFastaQual($$$){
  my($acc,$sff,$settings)=@_;
  my $basename=$acc;
  $basename=~s/_acc$//;

  my $command="sfffile -i '$acc' -o '$basename.sff' $sff";
  logmsg "COMMAND $command";
  system($command);
  die if $?;
  my $outlier_fastaqual = sff2fastaqual(["$$settings{tempdir}/newbler_outlier.sff"], $settings);
  
  return $outlier_fastaqual;
}

# creates qual and sequence fasta files for an SFF file (basecalling)
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
  #return $run_name; #debugging;
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
  #logmsg "Skipping newbler assembly for testing purposes";  return $run_name; #debugging
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
  
  #TODO remove all AMOS files before starting, so that they can be overwritten. Individual scripts have the force option but I cannot find the option in AMOScmp itself
	my $invoke_string = "AMOScmp";
	$invoke_string .= " -D TGT='$_'" for @afg_files;
	$invoke_string .= " -D REF='$_'" for @$ref_files;
	$invoke_string .= " $run_name";
	logmsg "Running $invoke_string";
	system($invoke_string); die if $?;

  # get CTG, AFG, and ACE file output too
  logmsg "Generating CTG, AFG, and ACE files.";
  $invoke_string="bank-report -b $$settings{tempdir}/amos_mapping.bnk CTG > $$settings{tempdir}/amos_mapping.ctg";
  logmsg $invoke_string;
  system($invoke_string);
  $invoke_string="amos2ace -o $$settings{tempdir}/amos_mapping.ace ";
  $invoke_string.="$_ " for @afg_files;
  $invoke_string.="$$settings{tempdir}/amos_mapping.ctg";
  logmsg $invoke_string;
  system($invoke_string);
  $invoke_string="toAmos -ace $$settings{tempdir}/amos_mapping.ace -o $$settings{tempdir}/amos_mapping.afg ";
  logmsg $invoke_string;
  system($invoke_string);

  #TODO sed 's/any trailing whitespace at the end of a line//' amos_mapping.ace

	return $run_name;
}

# combine newbler and velvet de novo assemblies
sub combineNewblerDeNovo($$$) {
	my ($newbler_basename, $settings) = @_;
	logmsg "Running Newbler combining on $newbler_basename";

  # get accessions for the reads that won't be in the final assembly: Outlier Repeat Singleton TooShort
	my (%outlier_reads, %repeat_reads,%singleton_reads,%tooShort_reads,$command);
	open(IN, '<', "$newbler_basename/assembly/454ReadStatus.txt") or die("Could not open 454ReadStatus.txt file\n");
	while (<IN>) {
		chomp;
		my ($read_id, $status) = split /\s+/;
		$outlier_reads{$read_id} = 1 if $status eq 'Outlier'; # problematic reads (chimeras, contamination, etc)
		$repeat_reads{$read_id} = 1 if $status eq 'Repeat'; # probably from repeat region, thus excluded from assembly
		$singleton_reads{$read_id} = 1 if $status eq 'Singleton'; # no overlaps found
		$tooShort_reads{$read_id} = 1 if $status eq 'TooShort'; # shorter than 50bp, or <15bp with paired end reads
	}
	close IN;

  # print out the accessions to file
	open(OUT, '>', "$$settings{tempdir}/newbler_outlier_acc") or die;
	print OUT "$_\n" for keys %outlier_reads;
	close OUT;
	open(OUT, '>', "$$settings{tempdir}/newbler_repeat_acc") or die;
	print OUT "$_\n" for keys %repeat_reads;
	close OUT;
	open(OUT, '>', "$$settings{tempdir}/newbler_singleton_acc") or die;
	print OUT "$_\n" for keys %singleton_reads;
	close OUT;
	open(OUT, '>', "$$settings{tempdir}/newbler_tooShort_acc") or die;
	print OUT "$_\n" for keys %tooShort_reads;
	close OUT;

  # create SFFs of the accessions not included in the assembly
  my($outlier_fastaqual,$repeat_fastaqual,$singleton_fastaqual,$tooShort_fastaqual);
  if(glob("$newbler_basename/sff/*.sff")){
    $outlier_fastaqual=sffFastaQual("$$settings{tempdir}/newbler_outlier_acc","$newbler_basename/sff/*.sff",$settings);
    $repeat_fastaqual=sffFastaQual("$$settings{tempdir}/newbler_repeat_acc","$newbler_basename/sff/*.sff",$settings);
    $singleton_fastaqual=sffFastaQual("$$settings{tempdir}/newbler_singleton_acc","$newbler_basename/sff/*.sff",$settings);
    $tooShort_fastaqual=sffFastaQual("$$settings{tempdir}/newbler_tooShort_acc","$newbler_basename/sff/*.sff",$settings);
  }
  # TODO if the input is a fasta file, then the SFF will not be present. Extract the relevant reads into fasta/qual files
  else{
    logmsg "Warning: SFF was not an input file.  Outlier, repeat, singleton, and tooShort reads will not be shown in their own files.";
    system("touch $$settings{tempdir}/newbler_repeat.sff.fasta $$settings{tempdir}/newbler_repeat.sff.qual $$settings{tempdir}/newbler_singleton.sff.fasta $$settings{tempdir}/newbler_singleton.sff.qual $$settings{tempdir}/newbler_tooShort.sff.fasta $$settings{tempdir}/newbler_tooShort.sff.qual $$settings{tempdir}/newbler_tooShort.sff.fasta $$settings{tempdir}/newbler_tooShort.sff.qual");
    die if $?;
  }

  # choose either Newbler contigs or scaffolds if they exist
  my $newblerAssembly="$newbler_basename/assembly/454AllContigs.fna";
  $newblerAssembly="$newbler_basename/assembly/454Scaffolds.fna" if(-s "$newbler_basename/assembly/454Scaffolds.fna" > 0);
  # begin the combining
	my $newbler_contigs = count_contigs($newblerAssembly); 
	my $combined_fasta_file = "$$settings{tempdir}/combined_in.fasta";
	my $numcontigs=0;
  system("cat '$newblerAssembly' '$$singleton_fastaqual[0]->[0]' > $combined_fasta_file");
  die if $?;
  system("toAmos -s '$combined_fasta_file' -o '$$settings{tempdir}/minimus.combined.afg'");
  system("minimus2 -D REFCOUNT=$numcontigs '$$settings{tempdir}/minimus.combined'");

	return "$$settings{tempdir}/minimus.combined.fasta";
}
# combine reference
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

	my $newbler_contigs = count_contigs("$newbler_basename/mapping/454AllContigs.fna");
	my $amos_contigs = count_contigs("$amos_basename.fasta");
	my $combined_fasta_file = "$$settings{tempdir}/combined_in.fasta";
	my $numcontigs=0;
	if($newbler_contigs < $amos_contigs){
		system("cat '$newbler_basename/mapping/454AllContigs.fna' '$amos_basename.fasta' '$$unmapped_fastaqual[0]->[0]' > $combined_fasta_file");
		$numcontigs=$newbler_contigs;
		logmsg("Newbler reference assembly selected as reference for minimus2");
	}
	else{
		system("cat '$amos_basename.fasta' '$newbler_basename/mapping/454AllContigs.fna' '$$unmapped_fastaqual[0]->[0]' > $combined_fasta_file");
		$numcontigs=$amos_contigs;
		logmsg("AMOScmp reference assembly selected as reference for minimus2");
	}
	die if $?;
	system("toAmos -s '$combined_fasta_file' -o '$$settings{tempdir}/minimus.combined.afg'");
	system("minimus2 -D REFCOUNT=$numcontigs '$$settings{tempdir}/minimus.combined'");

	return "$$settings{tempdir}/minimus.combined.fasta";
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
