#!/usr/bin/env perl

# run-assembly-pacbio: Perform standard assembly protocol operations on pacbio files
# Author: Lee Katz <lkatz@cdc.gov>

package PipelineRunner;
my ($VERSION) = ('$Id: $' =~ /,v\s+(\d+\S+)/o);

my $settings = {
  appname => 'cgpipeline',
  defaultGenomeSize=>5000000, # default genome: 5 MB. Many bacterial genomes fall around 5 MB (but there are many exceptions!)
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
use POSIX;

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH = $FSFind::LOG || *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

exit(main());

sub main() {
  $settings = AKUtils::loadConfig($settings);
  $$settings{assembly_min_contig_length} ||= 1000; #TODO put this into config file

  die(usage()) if @ARGV < 1;

  my @cmd_options = qw(keep tempdir=s outfile=s illuminaReads=s@ sffReads=s@ ccsReads=s@ expectedGenomeSize=i pacbiospec=s caspec=s);
  GetOptions($settings, @cmd_options) or die;
  $$settings{numcpus}||=AKUtils::getNumCPUs();
  if(!$$settings{expectedGenomeSize}){
    warn("WARNING: expected genome size was not given. I will set it to $$settings{defaultGenomeSize}");
    $$settings{expectedGenomeSize}=$$settings{defaultGenomeSize};
  }
  $$settings{ccsReads}||=[];

  my $stepNo=0;
  logmsg "STEP ".++$stepNo." receive input files";
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

  logmsg "STEP ".++$stepNo." convert H5 files to fastq";
  convertH5ToFastq(\@input_files,$settings);
  
  logmsg "STEP ".++$stepNo." modify input files and create FRG files";
  my $inputDir=highQualityReadsToFrg(\@input_files,$settings);

  logmsg "STEP ".++$stepNo." create FRG file for long reads";
  padLongReads($inputDir,$settings);

  logmsg "STEP ".++$stepNo." assembly";
  my $caAssembly=caAssembly($inputDir,$settings);

  logmsg "TODO STEP ".++$stepNo." AHA";
  my $ahaAssembly=ahaScaffolding(\@input_files,$caAssembly,$settings);

  logmsg "TODO STEP ".++$stepNo." dealing with unmapped reads";
  # STEP 6 unmapped short reads => de novo? or try to map again?

  my $final_seqs = AKUtils::readMfa($caAssembly);

  foreach my $seq (keys %$final_seqs) {
    delete $$final_seqs{$seq} if length($$final_seqs{$seq}) < $$settings{assembly_min_contig_length};
  }
  AKUtils::printSeqsToFile($final_seqs, $$settings{outfile}, {order_seqs_by_name => 1});

  logmsg "Output is in $$settings{outfile}";

  return 0;
}

########################
## Subs for major stages
########################

# Converts H5 to useful fastq files and adds the fq files
# to the ccsReads and @input_files arrays by reference.
sub convertH5ToFastq{
  my($input_files,$settings)=@_;
  my $newInputDir="$$settings{tempdir}/input";
  mkdir($newInputDir);
  my $longreadsInputDir="$newInputDir/longreads";
  mkdir("$longreadsInputDir");
  my $workingDir=AKUtils::mktempdir();
  return $newInputDir if(-f "$newInputDir/longreads.fastq" && -s "$newInputDir/longreads.fastq" > 100);

  my $convertQ=Thread::Queue->new(@$input_files);
  my @thr;
  for(0..ceil($$settings{numcpus}/2)){ # two threads per h5 file
    $thr[$_]=threads->new(sub{
      my($Q,$settings)=@_;
      my $bash5tools=AKUtils::fullPathToExec("bash5tools.py");
      my $i=0;
      while(defined(my $h5=$convertQ->dequeue)){
        my($b,$dir,$ext)=fileparse($h5,qw(.bas.h5 .h5));
        next if($ext!~/\.h5/i);
        my @c=("$bash5tools -f '$h5' -t CCS -s fastq -l 100 -q 0 -o '$workingDir/$b.ccs'",
               "$bash5tools -f '$h5' -t Raw -s fastq -l 2000 -q 0 -o '$workingDir/$b.raw'");
        command(\@c);
        command(["mv '$workingDir/$b.ccs.fastq' $newInputDir/",
          "mv '$workingDir/$b.raw.fastq' $longreadsInputDir/"]);
      }
    },$convertQ,$settings);
  }
  $convertQ->enqueue(undef) for(@thr);
  $_->join for(@thr);

  # grab all these new input files
  push(@{ $$settings{ccsReads} },glob("$newInputDir/*.ccs.fastq"));
  push(@$input_files,glob("$longreadsInputDir/*.raw.fastq"));
  @$input_files=grep(!/\.h5$/i,@$input_files); # remove h5 files

  return $newInputDir;
}

sub highQualityReadsToFrg{
  my($input_files,$settings)=@_;
  my $newInputDir="$$settings{tempdir}/input";
  mkdir($newInputDir);
  # STEP 2a standardize pacbio fastq files
  logmsg "Reading, standardizing, and cleaning long reads";
  my $longreadsFile="$newInputDir/longreads.fastq";
  if(!-e $longreadsFile){
    my $longreadsInputDir="$newInputDir/longreads";
    mkdir("$longreadsInputDir");
    for(my $i=0;$i<@$input_files;$i++){
      my $file=$$input_files[$i];
      my $filename=basename($file);
      my $standardizedFile="$longreadsInputDir/$filename";
      command("run_assembly_convertMultiFastqToStandard.pl -i '$file' -o '$standardizedFile' -m 200 2>&1") if(!-e $standardizedFile);
      $$input_files[$i]=$standardizedFile;
    }

    # filter and concatenate long reads at the same time
    my $command="run_assembly_trimClean.pl --min_avg_quality 6 --notrim --min_length 200 ";
    $command.="-i '$_' " for(@$input_files);
    $command.="-o $longreadsFile 2>&1 ";
    command($command);
    
    # One more round of filtering: keep longest reads only
    # This subroutine will edit the file in-place
    my $minRawReadLength=filterForCoverage($longreadsFile,$settings);
    command("mv $longreadsFile $longreadsFile.tmp.fastq");
    command("run_assembly_convertMultiFastqToStandard.pl -i $longreadsFile.tmp.fastq -o $longreadsFile -m $minRawReadLength");
    command("rm $longreadsFile.tmp.fastq");

    # Remove these temporary intermediate files.  They are huge and not necessary, even in a temp directory
    command("rm -rfv $longreadsInputDir");
  }

  # STEP 2b remove CCS reads from long reads (ie remove duplication)
  if(@{$$settings{ccsReads}}){
    logmsg "Reading and cleaning CCS pacbio reads";
    my $ccsreadsFile="$newInputDir/ccsreads.fastq";
    if(!-e $ccsreadsFile){
      my $command="run_assembly_trimClean.pl --min_avg_quality 20 --notrim --min_length 100 ";
      $command.="-i '$_' " for(@{$$settings{ccsReads}});
      $command.="-o $ccsreadsFile 2>&1";
      command($command);
    }
    @{ $$settings{ccsReads} }=($ccsreadsFile); # this command leaves us with one comprehensive ccsReads file
    removeReadsinCcsFromLongreads($$settings{ccsReads},$longreadsFile,$settings);
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
  # 454: sffToCA
  for my $sff (@{$$settings{sff}}){
    die("ERROR: sffToCA has not been implemented and so I cannot use the input file $sff");
    # suggested command: sffToCA -libraryname sff -output ${genome}454 $filename
    # where genome is the genome name and 454 indicates that this is the 454 platform, and filename is the path
  }

  return $newInputDir;
}

sub padLongReads{
  my($inputDir,$settings)=@_;
  my $longreadsFile="$inputDir/longreads.fastq";
  my $libraryname="longreadsCorrected";
  my $frg="$inputDir/$libraryname.frg";

  # STEP 3a create or locate appropriate spec files
  my $specFile=pacbiotocaSpecFile($inputDir,$settings);

  # STEP 3b pad the long reads with pacBioToCA
  #local $$settings{numcpus}=1; logmsg "DEBUGGING with numcpus=1";
  warn("NOTE: if there is any problem with pacBioToCA, please look at their development page at http://sourceforge.net/apps/mediawiki/wgs-assembler/index.php?title=PacBioToCA#Known_Issues\n");
  my $exec=AKUtils::fullPathToExec("pacBioToCA");
  # TODO remove overlap.sh if it exists, to help pacBioToCA continue if it was killed in the middle of running
  # TODO remove anything else to ensure continuity?
  command("rm '$inputDir/longreads.frg'",{warn_on_error=>1});
  my $command="(cd $inputDir; $exec -length 500 -partitions 200 -l $libraryname -t $$settings{numcpus} -noclean -s $specFile -fastq $longreadsFile $inputDir/*.frg)  2>&1 ";
  if(!-e $frg || -s $frg<10000){
    my $is_error=command($command);
  #  my $is_error=command($command,{warn_on_error=>1});
  #  if($is_error){
  #    logmsg "Trying pacBioToCA one more time with 1 cpu instead, which fixes a particular bug in some installations.";
  #    $command=~s/\s+\-t\s+\d+\s+/ -t 1 /;
  #    command($command);
  #  }
  } else { logmsg "Padded long reads have been finished already in $frg. Moving on.";}

  return $frg;
}
sub caAssembly{
  my($inputDir,$settings)=@_;
  my $finalAssembly="$$settings{tempdir}/assembly.contigs.fasta";
  # STEP 4a create or locate appropriate spec files
  my $specFile=runcaSpecFile($inputDir,$settings);
  # STEP 4b runCA
  my $caPrefix="$$settings{tempdir}/asm";
  my $asmPrefix="asm";

  # runCA is choking on abs paths
  $caPrefix=File::Spec->abs2rel($caPrefix);
  $specFile=File::Spec->abs2rel($specFile);
  $inputDir=File::Spec->abs2rel($inputDir);
  
  my $exec=AKUtils::fullPathToExec("runCA");

  #command("$exec -p $asmPrefix -d $caPrefix -s $specFile $inputDir/*.frg 2>&1"); # to include all reads
  command("$exec -p $asmPrefix -d $caPrefix -s $specFile $inputDir/longreads.frg 2>&1");
  # STEP 4c find the best assembly
  logmsg "Finished with the assembly! Grabbing your assembly now.";
  command("cp -v $caPrefix/9-terminator/asm.scf.fasta '$finalAssembly'");
  command("cp -v $caPrefix/9-terminator/asm.scf.qual '$finalAssembly.qual'");

  return $finalAssembly;
}

sub ahaScaffolding{
  my($input_files,$caAssembly,$settings)=@_;
  my @h5=grep(/\.h5/i,@$input_files);
  logmsg "AHA scaffolding is not complete, but if it were, then it could have been done using your H5 files and your assembly: ".join(", ",@h5,$caAssembly);
  return "";
}

####################
## Nitty-gritty subs
####################


# Find the coverage that an X bp min will give, and keep trying shorter reads until it reaches 20-40x.
# Arguments for settings: 
#  targetCoverage=>40 # the max coverage. The first level of coverage that is below this threshold will be used.
sub filterForCoverage{
  my($infile,$settings)=@_;
  my $targetCoverage=$$settings{targetCoverage}||50; # a maximum coverage

  # make a histogram of number of reads per length
  my %readLengthHist;
  open(IN,$infile) or die "Error: could not open $infile for reading: $!";
  while(<IN>){
    my $read=<IN>; chomp($read);
    my $readLength=length($read);
    $readLength=int(($readLength+449)/1000)*1000;
    $readLengthHist{$readLength}++;

    # burn two lines
    my $plusSign=<IN>; 
    my $qual=<IN>;
  }
  close IN;

  logmsg "Calculating minimum read length to coverage depth histogram";
  my @readLength=sort({$a<=>$b} keys(%readLengthHist));
  my %coverageHist;
  for(my $i=0;$i<@readLength;$i++){
    my $minReadLength=$readLength[$i];
    my $totalLength=0;
    for(my $j=$i;$j<@readLength;$j++){
      my $currReadLength=$readLength[$j];
      $totalLength+=($currReadLength*$readLengthHist{$currReadLength});
    }
    my $coverage=$totalLength/$$settings{expectedGenomeSize};
    $coverageHist{$minReadLength}=$coverage;
    print "$minReadLength\t$coverage\n";
  }
  
  # which coverage gives about the right target coverage?
  for my $readLength(@readLength){
    if($coverageHist{$readLength}<$targetCoverage){
      logmsg "$readLength yields a genome coverage of $coverageHist{$readLength} which is less than $targetCoverage. I'll go with that.";
      logmsg "Raw longreads coverage will be $coverageHist{$readLength}";
      return($readLength,$coverageHist{$readLength}) if wantarray;
      return $readLength;
    }
  }
  die "ERROR: Could not determine the right minimum read length cutoff";
}

sub ccsToFrg{
  my($fastq,$newInputDir,$settings)=@_;
  my $libraryname=basename($fastq,qw(.fastq .fq));
  my $frgFile="$newInputDir/$libraryname.frg";
  if(-e $frgFile && -s $frgFile>100){
    warn "NOTE: frg file $frgFile already exists. I will not overwrite it";
    return $frgFile;
  }

  # fastq to fasta/qual
  my ($fasta,$qual)=("$newInputDir/$libraryname.fna","$newInputDir/$libraryname.qual");
  logmsg "Creating a fasta/qual for $fastq";
  if(!-e $fasta || -s $fasta < 100){
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
  }

  # fasta/qual to frg
  my $command="convert-fasta-to-v2.pl -pacbio -s '$fasta' -q '$qual' -l $libraryname > '$frgFile'";
  command($command);

  return $frgFile;
}
sub illuminaToFrg{
  my($fastq,$newInputDir,$settings)=@_;
  my $libraryname=basename($fastq,qw(.fastq .fq));
  my $frgFile="$newInputDir/$libraryname.frg";
  if(-e $frgFile && -s $frgFile>100){
    warn "NOTE: frg file $frgFile already exists. I will not overwrite it";
    return $frgFile;
  }

  my $is_PE=AKUtils::is_fastqPE($fastq,$settings);

  # quality trim the illumina file
  my $cleanedFastq="$newInputDir/$libraryname.cleaned.fastq";
  # note: I removed the PE parameter because trimClean already detects if it is paired end
  command("run_assembly_trimClean.pl -i '$fastq' -o '$cleanedFastq' --min_quality 35 --min_avg_quality 30 --min_length 36 --bases_to_trim 20") if(!-e $cleanedFastq || -s $cleanedFastq <100);

  # make the frg
  my $command="fastqToCA -insertsize 300 20 -libraryname $libraryname -technology illumina -innie ";
  if($is_PE){
    $command.="-mates '$cleanedFastq' ";
  } else {
    $command.="-reads '$cleanedFastq' ";
  }
  $command.="1> $frgFile";
  command($command);
  return $frgFile;
}

sub pacbiotocaSpecFile{
  my ($inputDir,$settings)=@_;
  if($$settings{pacbiospec}){
    die "pacbio spec file not found!" if(!-f $$settings{pacbiospec});
    return $$settings{pacbiospec};
  }
  my $specFile="$inputDir/pacbiotoca.spec";
  return $specFile if(-e $specFile);

  command("cp '$FindBin::RealBin/../conf/pacbiotoca.spec' '$specFile' 2>&1");
  return $specFile;
}
sub runcaSpecFile{
  my ($inputDir,$settings)=@_;
  if($$settings{caspec}){
    die "CA spec file not found!" if(!-f $$settings{caspec});
    return $$settings{caspec};
  }
  my $specFile="$inputDir/runca.spec";
  return $specFile if(-e $specFile);

  command("cp '$FindBin::RealBin/../conf/runca.spec' '$specFile' 2>&1");
  return $specFile;
}

sub removeReadsinCcsFromLongreads{
  my ($ccsReadsArr,$longreadsFile,$settings)=@_;
  logmsg "TODO remove any reads found in the CCS reads from the long uncorrected reads file"
}

##########
## Utility
##########

# run a command
# settings params: warn_on_error (boolean)
# returns error code
# TODO run each given command in a new thread
sub command{
  my ($command,$settings)=@_;
  my $refType=ref($command);
  if($refType eq "ARRAY"){
    my $err=0;
    for my $c(@$command){
      $err+=command($c,$settings);
    }
    return $err;
  }
  logmsg "RUNNING COMMAND\n  $command";
  system($command);
  if($?){
    my $msg="ERROR running command $command\n  With error code $?. Reason: $!\n  in subroutine ".(caller(1))[3];
    if($$settings{warn_on_error}){
      logmsg $msg;
    } else {
      die $msg;
    }
  }
  return $?;
}

sub usage{
  "Usage: $0 input.fastq [, input2.fastq, ... inputX.h5 ...] [-o outfile.fasta] -c pacbio.ccs.fastq -i illumina.fastq -s 454.sff
  Input files should be the long filtered subreads from pacbio or raw h5 files
  -c, -i, -s
    You can use multiple files. 
    Designate each new file with a new -c, -i, or -s flag.
    -c is circular consensus sequence; -i is illumina; -s is SFF from 454 or iontorrent
    Multiple files from any platform are allowed.
    Paired end illumina data can be inputted using shuffled sequences in a singled file.

  -e expected genome size in bp
    This option is used for finding the best assembly after runCA finishes.
    Default: $$settings{defaultGenomeSize} but it is STRONGLY ENCOURAGED to supply a more exact size.
  -t tempdir
    Designate a different temporary directory
  -k
    If using default tempdir, use -k to keep temporary files
  --pacbiospec, --caspec to specify a different spec file for pacBioToCA or runCA
  "
}
