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
	$$settings{assembly_min_contig_length} ||= 500; #TODO put this into config file

	die("Usage: $0 input.sff [, input2.sff, ...] [-R references.mfa] [-C workdir]\n  Input files can also be fasta files.  *.fasta.qual files are considered as the relevant quality files") if @ARGV < 1;

	my @cmd_options = ('ChangeDir=s', 'Reference=s@', 'keep', 'tempdir=s', 'output=s');
	GetOptions($settings, @cmd_options) or die;

	$$settings{outfile} = $$settings{output};
	$$settings{outfile} ||= "$0.out.fasta";
	$$settings{outfile} = File::Spec->rel2abs($$settings{outfile});
  my($outfile,$outdir)=fileparse($$settings{outfile});
  system("mkdir -p $outdir");
	open(FH, '>', $$settings{outfile}) or die("Error writing to output file $$settings{outfile}");
	close FH;
	logmsg "Output file is " . $$settings{outfile} . "\n";

	my @ref_files = @{$$settings{Reference}} if defined $$settings{Reference};
	logmsg "No reference files supplied. Reverting to assembly mode" unless @ref_files;

	my @input_files = @ARGV;

	$$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($$settings{keep}));
	logmsg "Temporary directory is $$settings{tempdir}";
  system("mkdir -p $$settings{tempdir}");

	foreach my $file (@input_files, @ref_files) {
		$file = File::Spec->rel2abs($file);
		die("Input or reference file $file not found") unless -f $file;
	}

  my $fastaqualfiles=baseCall(\@input_files,$settings); #array of [0=>fna,1=>qual file] (array of arrays)

	my $final_seqs;

	if (@ref_files) {
		my $newbler_basename = runNewblerMapping(\@input_files, \@ref_files, $settings);

		my $amos_basename = runAMOSMapping($fastaqualfiles, \@ref_files, $settings);

		my $combined_filename = combineNewblerAMOS($newbler_basename, $amos_basename, $settings);

		$final_seqs = AKUtils::readMfa($combined_filename);
	} else {
    my($combined_filename);

		my $newbler_basename = runNewblerAssembly(\@input_files, $settings);
    my $newbler_assembly="$newbler_basename/assembly/454AllContigs.fna";
    $newbler_assembly="$newbler_basename/assembly/454Scaffolds.fna" if(-s "$newbler_basename/assembly/454Scaffolds.fna" > 0);

    $final_seqs=AKUtils::readMfa($newbler_assembly);
	}

	foreach my $seq (keys %$final_seqs) {
		delete $$final_seqs{$seq} if length($$final_seqs{$seq}) < $$settings{assembly_min_contig_length};
	}
	AKUtils::printSeqsToFile($final_seqs, $$settings{outfile}, {order_seqs_by_name => 1});

	logmsg "Output is in $$settings{outfile}";

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
  return 1 if(ext($file)=~/^(sff|fna|fasta|ffn|fa)$/i);
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
      die "The format of $file is incompatible. Acceptable formats are SFF and FASTA/QUAL";
    }
    push(@$fastaqualfiles,$seqQualPair);
  }
  return $fastaqualfiles;
}

# given a file containing one read accession per line and a pattern for all SFFs, extracts the reads and creates a new SFF
# the pattern for all SFFs is to be used in glob(), e.g. "sff/*.sff"
sub sffFastaQual($$$){
  my($acc,$sff,$settings)=@_;
  my $dir;
  ($acc,$dir)=fileparse($acc);
  my $basename=$acc;
  $basename=~s/_acc$//;

  my $command="sfffile -i '$$settings{tempdir}/$acc' -o '$$settings{tempdir}/$basename.sff' $sff";
  logmsg "COMMAND $command";
  system($command);
  die if $?;
  my $fastaqual = sff2fastaqual(["$$settings{tempdir}/$basename.sff"], $settings);
  
  return $fastaqual;
}

# creates qual and sequence fasta files for an SFF file (basecalling)
sub sff2fastaqual($$) {
	my ($sff_files, $settings) = @_;
	my @fastaqualfiles;
	foreach my $input_file (@$sff_files) {
		my ($sff_file, $sff_dir) = fileparse($input_file);
		my $invoke_string = "sffinfo -seq '$sff_dir/$sff_file' > '$$settings{tempdir}/$sff_file.fasta'";
		logmsg "Running $invoke_string";
		system($invoke_string) if(!-e "$$settings{tempdir}/$sff_file.fasta"); die "Failed when running\n  $invoke_string" if $?;
		$invoke_string = "sffinfo -qual '$sff_dir/$sff_file' > '$$settings{tempdir}/$sff_file.qual'";
		logmsg "Running $invoke_string";
		system($invoke_string) if(!-e "$$settings{tempdir}/$sff_file.qual"); die if $?;
		push(@fastaqualfiles, ["$$settings{tempdir}/$sff_file.fasta", "$$settings{tempdir}/$sff_file.qual"]);
	}
	return \@fastaqualfiles;
}
sub runNewblerMapping($$$) {
	my ($input_files, $ref_files, $settings) = @_;
	my $run_name = "$$settings{tempdir}/P__runMapping";
  #logmsg "Debugging"; return $run_name;
	logmsg "Executing Newbler mapping project $run_name";

	system("newMapping '$run_name'"); die if $?;
	foreach my $ref_file (@$ref_files) {
		system("setRef '$run_name' '$ref_file'"); die if $?;
	}
	foreach my $sff_file (@$input_files) {
		system("addRun '$run_name' '$sff_file'"); die if $?;
	}
	system("runProject '$run_name'"); die if $?;
  makeSpecialSffs($run_name,$settings);
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
	system("runProject '$run_name'"); die if $?;
  makeSpecialSffs($run_name,$settings);
	return $run_name;
}

sub runAMOSMapping($$$) {
	my ($input_files, $ref_files, $settings) = @_;
	my $run_name = "$$settings{tempdir}/amos_mapping";
  #logmsg "Debugging"; return $run_name;
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

sub makeSpecialSffs{
  my($newbler_basename,$settings)=@_;

  # make a hash of each read with its status
  my %reads;
  my $readStatus="$newbler_basename/assembly/454ReadStatus.txt";
  $readStatus="$newbler_basename/mapping/454ReadStatus.txt" if(!-e $readStatus);
  open(IN, '<', $readStatus) or die("Could not open $readStatus");
  while (<IN>) {
    chomp;
    my ($read_id, $status) = split /\s+/;
    $reads{$status}{$read_id}=1;
  }
  close IN;
  
  # print out the read IDs to files pertaining to their classification
  my @status=keys(%reads);
  for my $status (@status){
    my %specialReads=%{ $reads{$status} };
    open(OUT,'>',"$$settings{tempdir}/newbler_".$status."_acc") or die;
    print OUT "$_\n" for keys %specialReads;
    close OUT;
  }

  # create the SFF files
  if(glob("$newbler_basename/sff/*.sff")){
    for my $status (@status){
      sffFastaQual("$$settings{tempdir}/newbler_".$status."_acc","$newbler_basename/sff/*.sff",$settings);
    }
  } else {
    logmsg "Warning: SFF was not an input file.  Outlier, repeat, singleton, and tooShort reads will not be shown in their own files.";
    system("touch $$settings{tempdir}/newbler_repeat.sff.fasta $$settings{tempdir}/newbler_repeat.sff.qual $$settings{tempdir}/newbler_singleton.sff.fasta $$settings{tempdir}/newbler_singleton.sff.qual $$settings{tempdir}/newbler_tooShort.sff.fasta $$settings{tempdir}/newbler_tooShort.sff.qual $$settings{tempdir}/newbler_tooShort.sff.fasta $$settings{tempdir}/newbler_tooShort.sff.qual");
    die if $?;
  }

  return 1;
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

  my $combiningDir="$$settings{tempdir}/454combine";
  my $combined_filename="$combiningDir/minimus.combined.fasta";
  system("mkdir $combiningDir") if(!-d $combiningDir);
  die "Couldn't make combining directory because $!" if $?;
  system("run_assembly_combine.pl '$newbler_basename/mapping/454AllContigs.fna' '$amos_basename.fasta' -r '$$unmapped_fastaqual[0]->[0]' -t $combiningDir -o $combined_filename");
  die "Couldn't combine mapping assemblies" if $?;
  return $combined_filename;
}

