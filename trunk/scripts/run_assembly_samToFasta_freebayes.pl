#! /usr/bin/env perl

# Author: Lee Katz <lkatz@cdc.gov>
# Author: Eishita Tyagi <etyagi@cdc.gov>

package PipelineRunner;
my ($VERSION) = ('$Id: $' =~ /,v\s+(\d+\S+)/o);

my $settings = {
  appname  => 'cgpipeline',
  pileup_splitBases=>0,
  min_gap_size=>10,         # size of gap to split any contigs on. Fewer than this size is allowed to remain in a contig.
  pileup_min_frequency=>0.9,# percent of bases that must agree for a base call
};

use Data::Dumper;
use strict;
use warnings;
use Bio::Perl;
use FindBin;
use lib "$FindBin::RealBin/../lib";
$ENV{PATH} = "$FindBin::RealBin:".$ENV{PATH};
use AKUtils qw(logmsg);
use Getopt::Long;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use File::Spec;
use POSIX qw/floor ceil/;

use threads;
use Thread::Queue;

local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
exit(main());

sub main{
  $0=fileparse($0);
  my $settings=AKUtils::loadConfig($settings);
  die usage() if(@ARGV<1);

  GetOptions($settings,qw(assembly=s sam=s bam=s force tempdir=s keep outfile=s qual mask min_mapping_qual=i depth_allowedStdDeviations=i));
  $$settings{min_mapping_qual}||=0; # samtools mpileup -q parameter
  $$settings{numcpus}||=AKUtils::getNumCPUs();
  $$settings{depth_allowedStdDeviations}||=100;

  # check for required parameters
  for my $param (qw(assembly)){
    $$settings{$param} || die "Error: need $param parameter\n".usage();
  }
  # make some filenames absolute
  for my $param (qw(assembly sam bam)){
    next if(!$$settings{$param});
    $$settings{$param} = File::Spec->rel2abs($$settings{$param});
  }

  # set up the directory structure
  $$settings{tempdir}  ||= tempdir(File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($$settings{keep}));
  $$settings{tempdir}=File::Spec->rel2abs($$settings{tempdir});
  $$settings{outfile}||="$0.merged.fasta";
  logmsg "Temporary directory is $$settings{tempdir}";

  setupBuildDirectory($$settings{tempdir},$settings);

  my $bamIndex;
  my $bam=$$settings{bam};
    $bamIndex="$bam.bai" if($bam);
  my $sam=$$settings{sam};
  if($sam && !$bam){
    ($bam,$bamIndex)=createBam($sam,$settings);
  } else {
    my $samindexed="$$settings{assembly}.fai";
    my $alnSortedPrefix="$$settings{tempdir}/aln.sorted";
    my $alnSorted="$alnSortedPrefix.bam"; # bam file
    ($bam,$bamIndex)=indexBam($bam,$alnSortedPrefix,$alnSorted,$settings);
  }

  my $fastaOut=bamToFasta($bam,$bamIndex,$settings);
  #my ($newAssembly,$newAssemblyQual)=fastqToFastaQual($fastqOut,$settings);
  system("cp -v $fastaOut $$settings{outfile}");
  die "Could not cp $fastaOut to $$settings{outfile}:$!" if $?;
  logmsg "Output is in $$settings{outfile}";

  return 0;
}

sub setupBuildDirectory{
  my($tempdir,$settings)=@_;
  system("rm -rf $$settings{tempdir}") if($$settings{force});
  logmsg "Warning: temporary directory already exists. -f to force a fresh start" if(-e $$settings{tempdir});
  mkdir($$settings{tempdir});

  if($$settings{sam}){
    my $newSamPath="$tempdir/".basename($$settings{sam});
    system("ln -s $$settings{sam} $newSamPath") if(!-f $newSamPath);
    $$settings{sam}=$newSamPath;
  }

  my $newAssemblyPath="$tempdir/".basename($$settings{assembly});

  system("ln -s $$settings{assembly} $newAssemblyPath") if(!-f $newAssemblyPath);

  $$settings{assembly}=$newAssemblyPath;

  return 1;
}

# generate a .bam, .sorted.bam, .sorted.bam.bai and a database.fasta.fai file
sub createBam{
  my($sam,$settings)=@_;
  logmsg "Importing SAM to BAM";
  my $bamaln="$$settings{tempdir}/aln.bam";
  my $samindexed="$$settings{assembly}.fai";
  my $alnSortedPrefix="$$settings{tempdir}/aln.sorted";
  my $alnSorted="$alnSortedPrefix.bam"; # bam file
  my $bamIndex="$alnSorted.bai";

  # import the bam
  if(!-e $bamaln){
    logmsg "$bamaln does not exist. Creating it now.";
    command("samtools faidx $$settings{assembly}");
    command("samtools import $samindexed $sam $bamaln");
  } 

  # index the bam
  if(!-e $bamIndex){
    logmsg "Indexing the bam file";
    indexBam($bamaln,$alnSortedPrefix,$alnSorted,$settings);
  }

  return ($alnSorted,$bamIndex);
}

sub indexBam{
  my ($bamaln,$alnSortedPrefix,$alnSorted,$settings)=@_;
  my $bamIndex="$alnSorted.bai";
  return ($alnSorted,$bamIndex) if(-e $alnSorted && -e $bamIndex && -s $bamIndex > 100);
  command("samtools sort $bamaln $alnSortedPrefix") if(!-e $bamIndex || $bamIndex < 100);
  command("samtools index $alnSorted");
  return ($alnSorted,$bamIndex);
}

sub bamToFasta{
  my($bam,$bamIndex,$settings)=@_;

  my $mpileup="$$settings{tempdir}/mpileup.vcf";
  my $freebayes="$$settings{tempdir}/freebayes.vcf";
  my $basecallFile="$$settings{tempdir}/variants.unsorted.tsv";
  my $sortedBasecallFile="$$settings{tempdir}/variants.tsv";
  my $fastqOutStandard="$$settings{tempdir}/outStandard.fastq";
  my $fastqOut="$$settings{tempdir}/outSplitContigs.fastq";
  my($minDepth,$maxDepth)=covDepth($bam,$settings);

  if(!-e $freebayes){
    my $exec=AKUtils::fullPathToExec("freebayes");
    # output to a tmp file and then cp to the real destination, so that the file size doesn't matter when checking for the vcf
    command("$exec --left-align-indels --min-base-quality 20 --min-alternate-fraction 0.75 --min-coverage 5 -v '$freebayes.tmp' -b '$bam' -f '$$settings{assembly}'");
    command("cp '$freebayes.tmp' '$freebayes'");
  }
  
  # load the reference sequence
  my %consensus;
  my $ref=Bio::SeqIO->new(-file=>$$settings{assembly});
  while(my $seq=$ref->next_seq){
    $consensus{$seq->id}=[split(//,$seq->seq)];
  }
  $ref->close;

  # load up the alternate bases
  open(VCF,"<",$freebayes) or die "Could not open $freebayes: $!";
  while(<VCF>){
    next if(/^#/);
    chomp;
    my($rseq,$pos,$varId,$refBase,$altBase,$qual)=split /\t/;
    $consensus{$rseq}[$pos-1]=$altBase;
  }
  close VCF;

  # print the consensus to a file
  my $fastaOut="$$settings{tempdir}/consensus.fasta";
  open(FASTAOUT,">",$fastaOut) or die "Could not open $fastaOut for writing: $!";
  while(my($rseq,$bases)=each(%consensus)){
    $$bases[0]||="N"; # put at least some sequence here
    $_||="N" for(@$bases); # make sure that each base is declared
    my $consensus=join("",@$bases);
    $consensus=~s/(.{60})/$1\n/g;
    print FASTAOUT ">$rseq\n$consensus\n";
  }
  close FASTAOUT;
  return $fastaOut;
}

#########################
### utility methods
#########################

sub covDepth{
  my($bam,$settings)=@_;
  my $minimumAllowedCov=$$settings{minimumAllowedCov} || 5;
  my $allowedStdDeviations=$$settings{depth_allowedStdDeviations}; # how many stdevs to go out for depth?
  my $depthFile="$$settings{tempdir}/covDepth.tsv";
  logmsg "Finding the min/max depth within $allowedStdDeviations standard deviations";

  my @depth;
  if(!-f $depthFile || -s $depthFile < 1000){
    my $samtools=AKUtils::fullPathToExec("samtools");
    command("$samtools depth '$bam' > '$depthFile'");
  }
  
  open(DEPTH,$depthFile) or die "Could not open $depthFile for reading: $!";
  while(<DEPTH>){
    chomp;
    my ($rseq,$pos,$depth)=split /\t/;
    push(@depth,$depth);
  }
  close DEPTH;

  die "ERROR: There is no coverage depth information. This could happen if the wrong sam/bam and fasta files are paired together." if(!@depth);
  my $avgCov=average(\@depth);
  my $stdevCov=stdev(\@depth);

  my $max=ceil($avgCov+$allowedStdDeviations*$stdevCov);
  my $min=$avgCov-$allowedStdDeviations*$stdevCov;
  $min=$minimumAllowedCov if($min<$minimumAllowedCov);
  $min=floor($min);
  logmsg "The min/max is $min/$max, and the avg is $avgCov";

  return($min,$max);
}

sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}
sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}

# Prints lines as it comes out of a queue.  Enqueue undef to wrap it up.
# use this subroutine with a queue/threads
sub printer{
  my($file,$Q,$settings)=@_;
  logmsg "Preparing to write to $file";
  open(PRINTERFILE,">",$file) or die "Could not print to $file: $!"; 
  while(defined(my $line=$Q->dequeue)){
    print PRINTERFILE $line;
  }
  close PRINTERFILE;
  logmsg "Finished writing to $file";
  return 1;
}

sub command{
  my ($command)=@_;
  logmsg "RUNNING COMMAND\n  $command";
  system($command);
  die "ERROR running command $command\n  With error code $?" if($?);
  return 1;
}

sub usage{
  "This script converts a sam or bam to an assembly.
Usage: perl $0 -a reference.fasta -s assembly.sam -o assembly.fasta -q
  -s sam file
  -b bam file (don't use -b and -s in the same command)
  -a assembly reference file file
  -o (optional) final output file
    default: $0.merged.fasta
  -t (optional) This is where temporary files will be stored. 
    Default: /tmp/xxxxxx/ where xxxxxx is a random directory
  --min_mapping_qual An integer for the mapping quality needed.
    Default: 0
  -d allowed standard deviations of depth from the mean (3 std deviations to allow about 99% of all base calls, assuming a normal distribution)
    Default: 100 (REALLY permissive)

  No arguments should be given for the following options
  -f (optional) Force.
  -k (optional) keep temporary files around (about 2GB of files in my test run)
  -q to output quality files too. (assembly.fasta.qual)
  --mask to mask any bases that have a lower mapping quality or where the depth is not in 2 standard deviations of the mean
  ";
}

