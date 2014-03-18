#!/usr/bin/env perl

# run-assembly-illumina: Perform standard assembly protocol operations on illumina fastq file(s)
# Author: Lee Katz <lkatz@cdc.gov>
# Author: Eishita Tyagi <etyagi@cdc.gov>

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

	my @cmd_options = ('ChangeDir=s', 'Reference=s@', 'keep', 'tempdir=s', 'outfile=s', 'estimatedGenomeSize=i','numcpus=i', 'fast', 'concatenateWith=s');
	GetOptions($settings, @cmd_options) or die;
  $$settings{estimatedGenomeSize}||=5000000; # default: 5 MB
  $$settings{numcpus}||=1;

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
  mkdir($$settings{tempdir}) if(!-d $$settings{tempdir});
  $$settings{tempdir}=File::Spec->rel2abs($$settings{tempdir});
	logmsg "Temporary directory is $$settings{tempdir}";

	foreach my $file (@input_files, @ref_files) {
		$file = File::Spec->rel2abs($file);
		die("Input or reference file $file not found") unless -f $file;
	}

  my $fastqfiles=baseCall(\@input_files,$settings); #array of fastq files

	my ($combined_filename,$final_seqs);

	if (@ref_files) {
    my $smalt=AKUtils::fullPathToExec("smalt");
    my $samtools=AKUtils::fullPathToExec("samtools");
    my $bcfutils=AKUtils::fullPathToExec("bcftools");
    my $vcfutils=AKUtils::fullPathToExec("vcfutils.pl");

    logmsg "Concatenating the reference sequence(s) into one file and building a database";
    my $concatenatedReferences="$$settings{tempdir}/ref_seqs.fna";
    my $concatenatedReads="$$settings{tempdir}/reads.fastq";
    system("cat ".join(" ",@ref_files)." >$concatenatedReferences") if(!-e $concatenatedReferences); die "Could not concatenate references because $!" if $?;
    system("$smalt index -k 13 -s 2 $concatenatedReferences $concatenatedReferences") if(!-e "$concatenatedReferences.sma"); 
    die "Could not index the reference assembly $concatenatedReferences" if $?;

    # map each fastq against the database
    my @bam;
    for my $fastq (@$fastqfiles){
      my $sam="$$settings{tempdir}/".fileparse($fastq).".sam";
      my $bam="$$settings{tempdir}/".fileparse($fastq).".bam";
      my $sorted="$$settings{tempdir}/".fileparse($fastq).".sorted";

      my $isPe=AKUtils::is_fastqPE($fastq);
      logmsg "Mapping $fastq to assembly";
      if(!-e "$sorted.bam"){
        my $smaltxopts="-o $sam -n $$settings{numcpus}";
        if($isPe){
          my($mate1,$mate2)=deshuffleFastq($fastq,$settings);
          system("$smalt map $smaltxopts $concatenatedReferences $mate1 $mate2"); die "Error with smalt mapping" if $?;
        }else{
          system("$smalt map $smaltxopts $concatenatedReferences $fastq"); die "Error with smalt mapping" if $?;
        }
        system("$samtools view -bSh $sam > $bam"); die if $?;
        system("$samtools sort $bam $sorted"); die if $?;
        unlink($sam);unlink($bam);
      }
      push(@bam,"$sorted.bam");
    }
    logmsg "Done mapping";

    my $combinedBam="$$settings{tempdir}/combined.bam";
    if(@bam>1){
      logmsg "Combining output bam files.";
      system("$samtools merge $combinedBam ".join(" ",@bam)) if(!-e $combinedBam); die if $?;
    } else {
      system("cp $bam[0] $combinedBam"); die if $?;
    }

    # create the reference assembly from BAM
    logmsg "Changing mapped reads into an assembly";
    my $referenceAsm="$$settings{tempdir}/reference.$$.fasta";
    system("run_assembly_samToFasta.pl -b $combinedBam -a $concatenatedReferences -t $$settings{tempdir}/samToFasta -o $referenceAsm -q"); die if $?;

    # assemble unmapped reads
    # run de novo assembly, keeping unused reads
    logmsg "Recovering leftover unmapped reads";
    my $unmappedFastq="$$settings{tempdir}/unmapped";
    my ($uMate1,$uMate2)=($unmappedFastq."_1.fastq",$unmappedFastq."_2.fastq");
    system("bam2fastq -o '$unmappedFastq#.fastq' --no-aligned --unaligned --force $combinedBam"); die if $?;
    my $uShuffled=shuffleFastq([$uMate1,$uMate2],$settings);
    my $velvet_basename=runVelvetAssembly([$uShuffled],$settings);
    my $velvetAsm="$velvet_basename/contigs.fa";

    # mix in those de novo contigs
    $combined_filename="$$settings{tempdir}/referenceCombined.fasta";
    system("run_assembly_combine.pl -a $referenceAsm -a $velvetAsm -m 100 -o $combined_filename");
    # if combining doesn't work, just use the ref assembly
    if($?){
      logmsg "WARNING: Combining assemblies didn't work. Just using your reference assembly instead.";
      system("cp -v $referenceAsm $combined_filename");
    }
	} else {
    my($velvet_basename,$velvet_assembly,$spades_basename,$spades_assembly);
    $combined_filename="$$settings{tempdir}/denovoCombined.fasta";
    
    $spades_basename = runSpadesAssembly($fastqfiles,$settings);
    $spades_assembly = "$spades_basename/asm/contigs.fasta";

    $velvet_basename = runVelvetAssembly($fastqfiles,$settings);
    $velvet_assembly="$velvet_basename/contigs.fa";
    #system("run_assembly_combine.pl $spades_assembly $velvet_assembly -m 100 -o $combined_assembly -t $$settings{tempdir}/combine");

    # combine assemblies
    my $combined_assembly="$$settings{tempdir}/combine/combined_out.fasta";
    mkdir "$$settings{tempdir}/combine";
    mkdir "$$settings{tempdir}/combine/assembly";
    mkdir "$$settings{tempdir}/combine/reads";
    mkdir "$$settings{tempdir}/combine/ngsgam";
    system("ln -s $velvet_assembly $spades_assembly $$settings{tempdir}/combine/assembly/"); die if $?;
    system("ln -s ".join(" ".@$fastqfiles)." $$settings{tempdir}/combine/reads/"); die if $?;
    my $ngsgamCommand="run_assembly_combine_ngsgam.pl -asm $$settings{tempdir}/combine/assembly -reads $$settings{tempdir}/combine/reads -t $$settings{tempdir}/combine/ngsgam -e $$settings{expectedGenomeSize}";
    system($ngsgamCommand);
    if($?){
      logmsg "WARNING: either a problem with ngs-gam or it is not installed. Assemblies will not be merged.  Command was:\n  $ngsgamCommand" if $?;
      system("touch $combined_assembly"); die if $?;
    }

    # Even though the assemblies have been combined, find the best assembly.
    # Maybe the combined one isn't the best.
    system("run_assembly_chooseBest.pl $spades_assembly $velvet_assembly $combined_assembly -o $combined_filename");
    die "ERROR: run_assembly_chooseBest.pl failed" if $?;
	}

  $final_seqs = AKUtils::readMfa("$combined_filename");

	foreach my $seq (keys %$final_seqs) {
		delete $$final_seqs{$seq} if length($$final_seqs{$seq}) < $$settings{assembly_min_contig_length};
	}

  if(defined($$settings{concatenateWith})){
    my $sequence=join($$settings{concatenateWith},
                      values(%$final_seqs));
    $final_seqs={"concatenatedAssembly assembly_pipeline=CGP_$$settings{pipeline_version}" => $sequence};
  }
	AKUtils::printSeqsToFile($final_seqs, $$settings{outfile}, {order_seqs_by_name => 1});

	logmsg "Output is in $$settings{outfile}";

	return 0;
}

# determine file types
# determine file types
sub ext($){ # return the extension of a filename
  my($file)=@_;
  my($name,$path,$ext)=fileparse($file,qw(.fastq.gz .sff .fasta .fna .fa .fas .fastq .ffn .mfa));
  if(!$ext){
    die "Extension of $file could not be determined.";
  }
  $ext=~s/^.//; # remove that dot
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
  return 1 if(ext($file)=~/^(fastq|fq)$/i);
  return 0;
}
sub is_fastqGz($){
   my($file)=@_;
   return 1 if(ext($file)=~/^fastq\.gz$/i);
   return 0;
}
# wrapper for all base calling for input files
# TODO find out if there are mate pairs (observed-insert-lengths program via Velvet, or maybe by reading the headers?)
#TODO examine reads to see if the file overall belongs in "long" "short" etc: might just filter based on chemistry
sub baseCall($$){
  my($input_files,$settings)=@_;
  my $fastqfiles;
  foreach my $file (@$input_files){
    if(is_fastq($file)){
      push(@$fastqfiles,$file);
    }
    elsif(is_fastqGz($file)){
      my $targetFile=$$settings{tempdir}."/".fileparse($file);
      $targetFile=~s/\.gz$//;
      if(!-e $targetFile && -s $targetFile < 1){
        logmsg "have to decompress $file to $targetFile";
        system("gunzip -c '$file' > '$targetFile'");
      }
      push(@$fastqfiles,$targetFile);
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

sub runSpadesAssembly($$){
  my($fastqfiles,$settings)=@_;
  my $run_name = "$$settings{tempdir}/spades";
  
  # has it been run before?
  my $continueArg="";
  if(-d $run_name){
    $continueArg="--continue ";
  } else {
    system("mkdir -p $run_name");
    die if $?;
  }
  logmsg "Executing SPAdes assembly $run_name";
  
  my $spadesxopts=$$settings{spadesxopts}||"";
  my $fastqArg;
  my $carefulArg=""; # is set only if there are PE reads
  logmsg "WARNING: only using the first five input files for SPAdes" if(@$fastqfiles > 5);
  my @fastqfiles=("blurg",@$fastqfiles[0..4]); # unshift a bogus 0th element since we're in 1-based
     @fastqfiles=grep defined, @fastqfiles;    # remove extra elements if I grabbed too much
  for(my $i=1;$i<@fastqfiles;$i++){  # 1-based to match spades parameters
    if(AKUtils::is_fastqPE($fastqfiles[$i])){
      $fastqArg.="--pe$i-12 $fastqfiles[$i] ";
      $carefulArg="--careful";
    } else {
      $fastqArg.="-s $fastqfiles[$i] ";
    }
  }

  my $command="spades.py --tmp-dir $run_name/tmp -t $$settings{numcpus} -o $run_name/asm $fastqArg $carefulArg $continueArg";
  command($command);

  return $run_name;
}

sub runVelvetAssembly($$){
  my($fastqfiles,$settings)=@_;
  my $run_name = "$$settings{tempdir}/velvet";
  return $run_name if(-d $run_name);
  #system("mkdir -p $run_name") if(!-d $run_name); # is made by Voptimiser
  logmsg "Executing Velvet assembly $run_name";
  my $velvetPrefix=File::Spec->abs2rel("$$settings{tempdir}/auto"); # VelvetOptimiser chokes on an abs path

  # figure out how many threads will be used, even when considering memory
  # Estimate based on the first set of reads. This is just an estimate.
  # -109635 + 18977*ReadSize + 86326*GenomeSize + 233353*NumReads - 51092*K
  my $readLength=`head -6 $$fastqfiles[0] | tail -1 | wc -c`+0; # look at the second read
  my $numReads=`wc -l $$fastqfiles[0]`/4;
  my $memoryReqPerThread=-109635 + 18977*$readLength + 86326*$$settings{estimatedGenomeSize}/1000000 + 233353*$numReads/1000000 - 51092* sum(qw(19 21 23 25 27 29 31))/7;
  $memoryReqPerThread*=1024;
  my $totalMem=AKUtils::getFreeMem();
  my $memNumThreads=$totalMem/$memoryReqPerThread;
  my $numThreads=int(min($$settings{numcpus},$memNumThreads));
  if($numThreads<1){ # a zero will make it hang forever
    $numThreads=1;
    logmsg "WARNING you are low on memory. However, I will set the number of threads to 1 so that the script progresses.";
  }else{
    logmsg "I estimate that you have enough memory for $memNumThreads threads, and you have $$settings{numcpus} max threads. Threads set to $numThreads.";
  }

  my $step='-x 2'; #VelvetOptimiser hash step size. Cannot be odd or less than 2
  $step='-x 10' if($$settings{fast});

  my $command="VelvetOptimiser.pl -a -v -p $velvetPrefix -d $run_name -t $numThreads $step ";
  #warn "=======Debugging: only running hashes of 29 and 31\n"; $command.=" -s 29 -e 31 ";
  $command.="-f '";
  foreach my $fqFiles (@$fastqfiles){
    my $fileFormat="fastq"; # per the Velvet manual

    my $poly=AKUtils::is_fastqPE($fqFiles,{checkFirst=>20});
    my $readLength="short"; # per velvet docs
    if($poly==1){ # i.e. paired end
      $readLength="shortPaired";
    } elsif($poly>1){ die "Internal error with poly checker"; }

    $command.="-$readLength -$fileFormat $fqFiles ";
  }
  $command.="' 2>&1 ";
  logmsg "VELVET COMMAND\n  $command";
  system($command); die if $?;

  logmsg "Done. Velvet output is in $run_name";

  #TODO create dummy qual file (phred qual for each base is probably about 60-65). Or, if Velvet outputs it in the future, add in the qual parameter.

  # TODO make the amos2ace an option
  #logmsg "Creating ace file with amos2ace";
  #system("amos2ace $run_name/velvet_asm.afg"); logmsg "Warning: amos2ace failed" if $?; # make an ace file too

  return $run_name;
}

sub bowtie2{
  my($concatenatedReferences,$concatenatedReads,$settings)=@_;
  my $run_name="$$settings{tempdir}/bowtie2";
  mkdir $run_name if(!-d $run_name);

  my $sam="$run_name/bowtie2.sam";
  return $sam if(-e $sam && -s $sam > 10000);
  my $referenceLink="$run_name/ref.fasta";
  system("ln -sv $concatenatedReferences $referenceLink");
  die "ERROR Could not make a soft link for $concatenatedReferences to $referenceLink" if $?;
  system("bowtie2-build $referenceLink $referenceLink 2>&1");
  die "ERROR Could not create a bowtie2 index for $referenceLink" if $?;

  system("bowtie2 -p $$settings{num_cpus} -x $referenceLink -S $sam -q $concatenatedReads 2>&1");
  die "ERROR Could not run Bowtie2" if $?;

  return $sam;
}

# modified from Torsten Seemann <torsten@seemann.id.au>
# http://biostar.stackexchange.com/questions/14189/convert-fastq-to-afg/14198#14198
sub fastqToAfg{
  my($fastq,$outfile,$settings)=@_;
  $$settings{qOffset}||=33; # 33=sanger 64=illumina(<1.8)
  $$settings{aOffset}||=60;
  $$settings{afg_reportEvery}||=100000;

  return $outfile if(-s $outfile > 1000); # if the output file already has something in it, don't bother converting

  my $numLines=`wc -l $fastq`;
  my $numReads=$numLines/4;
  logmsg "Converting $numReads reads from $fastq to afg";

  open(FASTQ,$fastq) or die "Could not open $fastq because $!";
  open(AFG,">",$outfile) or die "Could not open $outfile for writing because $!";
  my $iid=0;
  while (my $eid = <FASTQ>) {

    die "bad fastq ID line '$eid'\n".<FASTQ> unless $eid =~ m/^\@/;
    $eid = substr $eid, 1;  # remove '@'

    my $seq = scalar(<FASTQ>);
    chomp $seq;
    $seq=~s/(.{60})/\1\n/g; #newline every 60 bp

    my $id2 = scalar(<FASTQ>);
    die "bad fastq ID2 line '$id2'" unless $id2 =~ m/^\+/;

    my $qlt = scalar(<FASTQ>);
    chomp $qlt;

    $iid++;
    print AFG "{RED\n";
    print AFG "iid:$iid\n";
    print AFG "eid:$eid";  # already has \n
    print AFG "seq:\n$seq\n.\n";
    print AFG "qlt:\n";
    my @q = split m//, $qlt;
    @q = map { chr( ord($_)-$$settings{qOffset} + $$settings{aOffset} ) } @q;
    my $qltAfg=join("",@q);
    $qltAfg=~s/(.{60})/\1\n/g; #newline every 60 bp
    print AFG $qltAfg."\n.\n";
    print AFG "}\n";
    if($iid % $$settings{afg_reportEvery} == 0){
      $|++;
      print ".";
      $|--;
    }
    print "." if($iid % $$settings{afg_reportEvery} == 0);
    if(0 && $iid>20){
      logmsg "DEBUG";
      last;
    }
  }
  print "\n";
  close FASTQ; close AFG;
  return $outfile;
}

# quick and dirty
sub fastqToFasta{
  my($fastq,$outPrefix,$settings)=@_;
  $$settings{qOffset}||=33;
  my $outfasta="$outPrefix.fasta";
  my $outqual="$outPrefix.qual";
  return ($outfasta,$outqual) if(-e $outfasta && -e $outqual);

  # make things fast by getting the length of the read ahead of time
  my $secondLine=`head -2 $fastq|tail -1`;
  chomp($secondLine);
  my $readLength=length($secondLine);
  die "The read length has been determined to be $readLength, which is too small" if($readLength<10);

  open(FASTQ,$fastq) or die "Could not open $fastq for reading because $!";
  open(FASTA,">",$outfasta) or die "Could not write to $outfasta because $!";
  open(QUAL,">",$outqual) or die "Could not write to $outqual because $!";
  my $i=0;
  while(<FASTQ>){
    my $whichLine=$i%4;
    if($whichLine==0){
      my $defline=$_;
      $defline=~s/^\@/>/;
      $defline=~s/^\s+|\s+$//g;
      print FASTA $defline."\n";
      print QUAL $defline."\n";
    } elsif ($whichLine==1){
      print FASTA $_;
    } elsif ($whichLine==3){
      my $qual=$_;
      chomp $qual;
      for(my $j=0;$j<$readLength;$j++){
        my $newQual=(ord(substr($qual,$j,1))-$$settings{qOffset})." ";
        print QUAL $newQual;
      }
      print QUAL "\n";
    }
    $i++;
  }
  close FASTQ; close FASTA; close QUAL;
  return ($outfasta,$outqual);
}

sub shuffleFastq{
  my($fastqs,$settings)=@_;
  my $fastq="$$settings{tempdir}/shuffled.$$.fastq";
  my($mate1,$mate2)=@$fastqs;
  open(MATE1,$mate1) or die "Could not open $mate1:$!";
  open(MATE2,$mate2) or die "Could not open $mate2:$!";
  open(SHUFFLED,">",$fastq) or die "Could not open $fastq:$!";
  while(my $id=<MATE1>){
    my $read1=$id;
    my $read2;
    $read1.=<MATE1> for(1..3);
    $read2.=<MATE2> for(1..4);
    print SHUFFLED "$read1$read2";
  }
  close SHUFFLED; close MATE1; close MATE2;
  return $fastq;
}
sub deshuffleFastq{
  my($fastq,$settings)=@_;

  # random names for the two mates
  my @set = ('0' ..'9', 'A' .. 'F');
  my $rand; $rand.=join("",map($set[rand @set],1..8));
  my $mate1="$$settings{tempdir}/$rand.$$.1.fastq";
  my $mate2="$$settings{tempdir}/$rand.$$.2.fastq";

  my $i=0;
  open(MATE1,">",$mate1) or die "Could not open mate1 for writing:$!";
  open(MATE2,">",$mate2) or die "Could not open mate2 for writing:$!";
  open(FASTQ,$fastq) or die "Could not open interleved fastq $fastq: $!";
  while(<FASTQ>){
    my $mod=$i%8;
    if($mod<4){ # 0,1,2,3 belong to mate1
      print MATE1 $_;
    } else {    # 4,5,6,7 belong to mate2
      print MATE2 $_;
    }
    $i++;
  }
  close FASTQ;close MATE1; close MATE2;
  return($mate1,$mate2);
}

sub bamUnmappedReads{
  my($bam,$settings)=@_;
  my $fastq="$$settings{tempdir}/unmapped.fastq";
  logmsg "Getting unmapped reads from $bam and putting them into $fastq";
  open(BAM,"samtools view $bam|") or die "Could not open $bam: $!";
  open(FASTQ,">",$fastq) or die "Could not open $fastq for writing:$!";
  while(<BAM>){
    my $flags=samFlags($_,$settings);
    next if(!$$flags{unmapped});
    chomp;
    my($id,$sequence,$qual)=(split(/\t/,$_))[0,9,10];
    print FASTQ "\@$id\n$sequence\n+\n$qual\n";
  }
  close FASTQ;close BAM;
  return $fastq;
}

# return which flags are true/false in a sam line
sub samFlags{
  my ($line,$settings)=@_;
  my $f=(split(/\t/,$line))[1];
  #my $dec=dec2bin($f);
  my %flag=(dup=>0,poorQual=>0,notPrimary=>0,mate2=>0,mate1=>0,mateReversed=>0,reversed=>0,mateUnmapped=>0,unmapped=>0,mappedWithMate=>0,paired=>0);
  if($f-0x0400 > -1){
    $flag{dup}=1;
    $f-=0x0400;
  }
  if($f-0x0200 > -1){
    $flag{poorQual}=1;
    $f-=0x0200;
  }
  if($f-0x0100 > -1){
    $flag{notPrimary}=1;
    $f-=0x0100;
  }
  if($f-0x0080 > -1){
    $flag{mate2}=1;
    $f-=0x0080;
  }
  if($f-0x0040 > -1){
    $flag{mate1}=1;
    $f-=0x0040;
  }
  if($f-0x0020 > -1){
    $flag{mateReversed}=1;
    $f-=0x0020;
  }
  if($f-0x0010 > -1){
    $flag{reversed}=1;
    $f-=0x0010;
  }
  if($f-0x0008 > -1){
    $flag{mateUnmapped}=1;
    $f-=0x0008;
  }
  if($f-0x0004 > -1){
    $flag{unmapped}=1;
    $f-=0x0004;
  }
  if($f-0x0002 > -1){
    $flag{mappedWithMate}=1;
    $f-=0x0002;
  }
  if($f-0x0001 > -1){
    $flag{paired}=1;
    $f-=0x0001;
  }
  if($f>0){
    die "ERROR: there is something unaccounted for in the flag for this line: $line\nFound flags:\n".Dumper(\%flag);
  }
  return \%flag;
}

# decimal to binary
sub dec2bin {
  my $str = unpack("B32", pack("N", shift));
  $str =~ s/^0+(?=\d)//;   # otherwise you'll get leading zeros
  return $str;
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
	"Usage: $0 input.fastq [, input2.fastq, ...] [-o outfile.fasta] [-R references.mfa] [-t tempdir]
  Input files can also be fasta files.  *.fasta.qual files are considered as the relevant quality files
  -concat LINKER to concatenate the contigs with a linker, e.g. NNNNNCACACACTTAATTAATTAAGTGTGTGNNNNN
  -n numcpus Default: 1
  -e estimated genome size, e.g., -e 3000000 for a 3MB genome
  "
}
