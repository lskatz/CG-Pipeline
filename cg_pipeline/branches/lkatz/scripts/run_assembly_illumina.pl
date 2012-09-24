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
#use Bio::Tools::Run::BWA;
#use Bio::Tools::Run::Samtools;

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH = $FSFind::LOG || *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

exit(main());

sub main() {
	$settings = AKUtils::loadConfig($settings);
	$$settings{assembly_min_contig_length} ||= 100; #TODO put this into config file

	die(usage()) if @ARGV < 1;

	my @cmd_options = ('ChangeDir=s', 'Reference=s@', 'keep', 'tempdir=s', 'outfile=s');
	GetOptions($settings, @cmd_options) or die;

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
    #die "Using references for illumina assembly is not working right now";
    logmsg "Warning: paired end reference illumina assemby is not supported right now (I'm not sure if you supplied a paired end parameter though)";

    logmsg "Concatenating the reference sequence(s) into one file. And then all reads into one file.";
    my $concatenatedReferences="$$settings{tempdir}/ref_seqs.fna";
    my $concatenatedReads="$$settings{tempdir}/reads.fastq";
    system("cat ".join(" ",@ref_files)." >$concatenatedReferences") if(!-e $concatenatedReferences); die "Could not concatenate references because $!" if $?;

    # create a 0-byte file for reads, and then append reads
    if(!-e $concatenatedReads){
      system("echo -n '' > $concatenatedReads"); die "Could not create file $concatenatedReads" if $?;
      for(@$fastqfiles){
        if(is_fastqGz($_)){
          logmsg "Gunzipping/concating $_";
          system("gunzip -c '$_' >> $concatenatedReads");
        } elsif(is_fastq($_)){
          logmsg "Concating $_";
          system("cat '$_' >> $concatenatedReads")
        }
        die "Could not concatenate reads because $!" if $?;
      }
    }
    
    # bowtie2
    # (for p in $genomes2; do echo "Mapping $p"; run_pipeline create -p $p; mkdir -p $p/build/assembly/2010EL-1786.fasta/bowtie2; ref="$p/build/assembly/2010EL-1786.fasta/bowtie2/ref.fa"; ln -sv `pwd -P`/../refGenomes/2010EL-1786.fasta $ref; bowtie2-build $ref $ref; gunzip -c combined/$p.filtered.fastq.gz > combined/$p.filtered.fastq; bowtie2 -p 12 -x $ref -S "$p/build/assembly/2010EL-1786.fasta/bowtie2/mapping.sam" -q combined/$p.filtered.fastq; rm combined/$p.filtered.fastq; run_assembly_samToFasta.pl -s "$p/build/assembly/2010EL-1786.fasta/bowtie2/mapping.sam" -a $ref -t "$p/build/assembly/2010EL-1786.fasta/bowtie2/toFasta" -o "$p/build/assembly/2010EL-1786.fasta/bowtie2/bowtie.fasta"; echo "Finished with $p"; done;) 2>&1 | tee --append bowtie2_again.log
    die "TODO";
    my $bowtie_sam=bowtie2($concatenatedReferences,$concatenatedReads,$settings);
    my $bowtie_assembly=samToFasta($bowtie_sam,$concatenatedReferences,$settings);
    my $bowtie_unmappedContigs=assembleUnmappedContigs($bowtie_sam,$settings);

    # should we concatenate unmapped contigs onto the reference assembly?  I think so, because contamination can be filtered out later.

    $combined_filename="$$settings{tempdir}/mappingAssembly_combined.fasta";

    
	} else {
    my($velvet_basename,$velvet_assembly);

    $velvet_basename = runVelvetAssembly($fastqfiles,$settings);
    $velvet_assembly="$velvet_basename/contigs.fa";

    # in case assemblies from multiple assemblers are combined
    $combined_filename=$velvet_assembly;
	}

  $final_seqs = AKUtils::readMfa("$combined_filename");

	foreach my $seq (keys %$final_seqs) {
		delete $$final_seqs{$seq} if length($$final_seqs{$seq}) < $$settings{assembly_min_contig_length};
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

sub runVelvetAssembly($$){
  my($fastqfiles,$settings)=@_;
  my $run_name = "$$settings{tempdir}/velvet";
  #system("mkdir -p $run_name") if(!-d $run_name); # is made by Voptimiser
  logmsg "Executing Velvet assembly $run_name";
  my $velvetPrefix=File::Spec->abs2rel("$$settings{tempdir}/auto"); # VelvetOptimiser chokes on an abs path
  my $command="VelvetOptimiser.pl -a -v -p $velvetPrefix -d $run_name ";
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
sub nesoni{
  my($concatenatedReferences,$concatenatedReads,$settings)=@_;
  logmsg "Running nesoni samshrimp";
  my $nesoni_run="$$settings{tempdir}/nesoni";
  return "$nesoni_run/contigs.fasta" if(-s "$nesoni_run/contigs.fasta" > 100);
  mkdir($nesoni_run) if(!-d $nesoni_run);

  system("nesoni samshrimp: $nesoni_run $concatenatedReferences reads: $concatenatedReads --sam-unaligned no") unless (-e "$nesoni_run/alignments.bam");
  if($?){
    logmsg "Problem with Nesoni and/or SHRiMP. Exiting nesoni.";
    system("touch $nesoni_run/contigs.fasta");
    return "$nesoni_run/contigs.fasta";
  }
  logmsg "Running nesoni samconsensus";
  system("nesoni samconsensus: $nesoni_run --majority 0.7 --cutoff 0.7") unless(-e "$nesoni_run/consensus.fa");
  die "Problem with Nesoni and/or Samtools" if $?;
  
  # read the consensus sequence into contigs
  my $unmaskedContigs="$nesoni_run/consensus.fa";
  my $contigs={};
  my $c=0;
  my $mappings=AKUtils::readMfa($unmaskedContigs,$settings);
  for my $contig(values(%$mappings)){
    my @contigs=split(/N+/,$contig);
    @contigs=grep(!/^\s*$/,@contigs);
    $$contigs{++$c}=$_ for(@contigs);
  }

  AKUtils::printSeqsToFile($contigs,"$nesoni_run/contigs.fasta");

  return "$nesoni_run/contigs.fasta";
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
  "
}
