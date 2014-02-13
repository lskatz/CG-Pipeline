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
  my $sortedAsm=sortAssemblies(\@asm,$settings); # best to worst
  my $gamDir=gamCreate($sortedAsm,$asmBam,$settings);
  my $fasta=gamMerge($gamDir,$sortedAsm,$settings);

  system("cat $fasta"); die if $?;

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
      #$sorted="$$settings{tempdir}/".fileparse($assembly).".sorted";
      system("samtools merge ".join(" ",@bam)." $mappingTempdir/merged.bam");
      die if $?;
      system("samtools sort $mappingTempdir/merged.bam $mappingTempdir/".fileparse($assembly));
      die if $?;
    } else {
    # just sort the one bam
      system("samtools sort $bam[0] $sorted");
      die if $?;
    }
    
    # also need to index the bam file
    system("samtools index $sorted.bam");
    die if $?;
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
  my($sortedAsm,$asmBam,$settings)=@_;

  logmsg "Running gam-create";

  my $gamDir="$$settings{tempdir}/gam";
  mkdir $gamDir if(!-d $gamDir);

  my %insRange;
  for my $asm(@$sortedAsm){
    my $bam=$$asmBam{$asm};
    # find out the insert size
    my $readMetrics=_readMetrics($bam,$settings);
    if($$readMetrics{medianFragmentLength} && $$readMetrics{medianFragmentLength}=~/(\d+)\[(\d+),(\d+)\]/){
      my $median=$1;
      my $q1Dist=abs($median-$2);
      my $q2Dist=abs($3-$median);
      my $max=$median + 3*$q2Dist;
      my $min=$median - 3*$q1Dist;
      $min=1 if($min<1);
      $insRange{$asm}=[$min,$max];
    } else {
      die "ERROR: Could not determine the fragment size from $bam";
    }
  }

  # create the master/slave files needed for GAM, with min/max insert sizes
  my $masterGam="$gamDir/master.PE.bams.txt";
  my $slaveGam ="$gamDir/slave.PE.bams.txt";
  # create the master file with the bam of the best assembly
  open(MASTER,">",$masterGam) or die "ERROR: Could not open $masterGam for writing: $!";
  print MASTER "$$asmBam{$$sortedAsm[0]}\n".join(" ",@{$insRange{$$sortedAsm[0]}})."\n";
  close MASTER;

  # slave file with the other assemblies
  open(SLAVE,">",$slaveGam) or die "ERROR: Could not open $slaveGam for writing: $!";
  for(my $i=1;$i<@$sortedAsm;$i++){
    print SLAVE "$$asmBam{$$sortedAsm[$i]}\n".join(" ",@{$insRange{$$sortedAsm[$i]}})."\n";
  }
  close SLAVE;

  system("gam-create --master-bam $masterGam --slave-bam $slaveGam  --min-block-size 10 --output $gamDir/out");
  die "ERROR: Problem with gam-create" if $?;
  return $gamDir;
}

sub gamMerge{
  my($gamDir,$sortedAsm,$settings)=@_;
  my $outPrefix="$gamDir/gam-ngs.out";
  logmsg "Running gam-merge";
  my $command="gam-merge --blocks-file $gamDir/out.blocks --master-fasta $$sortedAsm[0] --slave-fasta $$sortedAsm[1] --master-bam $gamDir/master.PE.bams.txt --slave-bam $gamDir/slave.PE.bams.txt --min-block-size 10 --threads $$settings{numcpus} --output $outPrefix 1>&2";
  logmsg $command;
  system($command);
  die "ERROR: problem with gam-merge" if $?;
  return "$outPrefix.gam.fasta";
}

#### UTIL subroutines
sub _readMetrics{
  my($reads,$settings)=@_;
  my $metrics=`run_assembly_readMetrics.pl --fast '$reads'`;
  die "ERROR: problem with run_assembly_readMetrics.pl: $!" if $?;
  chomp($metrics);
  my ($header,$values)=split(/\n/,$metrics);
  my @header=split /\t/,$header;
  my @values=split /\t/,$values;
  my %metrics;
  @metrics{@header}=@values;
  return \%metrics;
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
