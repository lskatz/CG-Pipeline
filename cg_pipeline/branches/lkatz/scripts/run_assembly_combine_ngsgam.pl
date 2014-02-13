#!/usr/bin/env perl

# Combines two or more assemblies with ngs-gam. Requires Illumina PE reads.

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use File::Temp qw/tempdir/;
use FindBin;
use lib "$FindBin::RealBin/../lib";
$ENV{PATH} = "$FindBin::RealBin:".$ENV{PATH};
use AKUtils;

$0=fileparse $0;
sub logmsg {print STDERR "$0: ".(caller(1))[3].": @_\n";}

exit(main());
sub main{
  my $settings={};
  GetOptions($settings,qw(help asmDir=s readsDir=s tempdir=s keep numcpus=i expectedGenomeSize=i)) or die $!;
  die usage() if($$settings{help});
  $$settings{numcpus}||=1;
  $$settings{expectedGenomeSize}||=0;
  $$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($$settings{keep}));
  mkdir $$settings{tempdir} if(!-d $$settings{tempdir});
  logmsg "Temporary directory is $$settings{tempdir}";
  
  # error check
  for(qw(asmDir readsDir)){
    die "ERROR: need $_\n".usage() if(!$$settings{$_});
    die "ERROR: $_ is not a directory\n".usage() if(!-d $$settings{$_});
  }

  # check for assemblies and reads
  my @asm=glob("$$settings{asmDir}/*.{fasta,fna,fa}");
  die "ERROR: no assemblies found in $$settings{asmDir} \n".usage() if(!@asm);
  my @read=(glob("$$settings{readsDir}/*.{fastq,fastq.gz}"));
  die "ERROR: no assemblies found in $$settings{readsDir} \n".usage() if(!@read);

  my $merged;
  my $asmBam=mapReads(\@read,\@asm,$settings);
  @asm=sortAssemblies(\@asm,$settings); # best to worst
  gamCreate();
  gamMerge();

  return 0;
}

sub mapReads{
  my($readArr,$asm,$settings)=@_;
  my @read=@$readArr;
  
  # map each read set against each assembly
  my $mappingTempdir="$$settings{tempdir}/mapping";
  mkdir $mappingTempdir;
  my %mergedBam;
  for my $assembly(@$asm){
    my $sorted="$$settings{tempdir}/".fileparse($assembly).".sorted";
    $mergedBam{$assembly}="$sorted.bam";
    if(-f "$sorted.bam" && -s "$sorted.bam" > 1){
      logmsg "Found $sorted.bam. Not recreating.";
      next;
    }

    system("smalt index $assembly $assembly");
    die if $?;
    # map all PE read sets
    my @bam;
    for my $reads(@read){
      next if(!AKUtils::is_fastqPE($reads));
      my $b="$mappingTempdir/".fileparse($reads);
      my $bam="$b.bam";
      my $reads1="$mappingTempdir/reads1.fastq";
      my $reads2="$mappingTempdir/reads2.fastq";

      system("run_assembly_shuffleReads.pl -d $reads 1> $reads1 2> $reads2");
      die "ERROR: problem with run_assembly_shuffleReads.pl" if $?;
      system("smalt map -n $$settings{numcpus} -f samsoft $assembly $reads1 $reads2 | samtools view -bS -T $assembly - > $bam");
      die if $?;
      push(@bam,$bam);
    }
    
    # merge and sort all reads' bams if there are >1
    if(@read > 1){
      $sorted="$$settings{tempdir}/".fileparse($assembly).".sorted";
      system("samtools merge ".join(" ",@bam)." $mappingTempdir/merged.bam");
      die if $?;
      system("samtools sort $mappingTempdir/merged.bam $mappingTempdir/".fileparse($assembly));
      die if $?;
    } else {
    # just copy the reads if there is only one bam
      system("cp -v $read[0] $sorted.bam");
      die if $?;
    }
  }
  return \%mergedBam;
}

# Sort the assemblies by assembly score
sub sortAssemblies{
  my($asm,$settings)=@_;
  my %asmScore;
  for(@$asm){
    my $command="run_assembly_metrics.pl -number -s assemblyScore $_ ";
    $command.="-e $$settings{expectedGenomeSize} " if($$settings{expectedGenomeSize});
    my $score=`$command` + 0;
    $asmScore{$_}=$score;
  }
  my @asm=sort{$asmScore{$b} <=> $asmScore{$a}} @$asm;

  return \@asm;
}

sub gamCreate{
  # gam-create --master-bam master.PE.bams.txt --slave-bam slave.PE.bams.txt --min-block-size 10
}

sub gamMerge{
  # gam-merge --blocks-file out.blocks --master-fasta spades.fasta --master-bam master.PE.bams.txt --slave-fasta velvet.fasta --slave-bam slave.PE.bams.txt --min-block-size 10 --threads 32 --output gam-ngs.out 2> gam-merge.err
}

sub usage{
  "Combines two or more assemblies with ngs-gam package.
  Usage: $0 -asm asmDir -reads readDir > out.fasta
    -asm asmDir A directory of assemblies
    -reads readsDir A directory of cleaned and shuffled PE reads
    -t tempDir
    -e expectedGenomeSize, e.g. for 3MB, -e 3000000
  "
}
