#! /usr/bin/env perl

# Original concept by Eishita Tyagi
# Author: Eishita Tyagi <etyagi@cdc.gov>
# Author: Lee Katz <lkatz@cdc.gov>

package PipelineRunner;
my ($VERSION) = ('$Id: $' =~ /,v\s+(\d+\S+)/o);

my $settings = {
  appname  => 'cgpipeline',
  pileup_splitBases=>0
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

exit(main());

sub main{
  $0=fileparse($0);
  my $settings=AKUtils::loadConfig($settings);
  die usage() if(@ARGV<1);

  GetOptions($settings,qw(assembly=s sam=s force tempdir=s keep outfile=s qual));

  # check for required parameters
  for my $param (qw(assembly sam)){
    $$settings{$param} || die "Error: need $param parameter\n".usage();
  }
  # make some filenames absolute
  for my $param (qw(assembly sam)){
    $$settings{$param} = File::Spec->rel2abs($$settings{$param});
    die "Error: $$settings{$param} doesn't exist!" unless(-e $$settings{$param});
  }

  # set up the directory structure
  $$settings{tempdir}  ||= tempdir(File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($$settings{keep}));
  $$settings{tempdir}=File::Spec->rel2abs($$settings{tempdir});
  $$settings{outfile}||="$0.merged.fasta";
  logmsg "Temporary directory is $$settings{tempdir}";

  my $sam=$$settings{sam};

  # start the computation

  setupBuildDirectory($$settings{tempdir},$settings);

  my ($bam,$bamIndex)=createBam($sam,$settings);
  my $fastqOut=bamToFastq($bam,$bamIndex,$settings);
  my ($newAssembly,$newAssemblyQual)=fastqToFastaQual($fastqOut,$settings);

  open(OUT,">$$settings{outfile}") or die "Could not open $$settings{outfile} because $!";
  print OUT $newAssembly;
  close OUT;
  open(OUT,">$$settings{outfile}.qual") or die "Could not open $$settings{outfile}.qual because $!";
  print OUT $newAssemblyQual;
  close OUT;
  logmsg "Output is in $$settings{outfile}";

  return 0;

  # print out the new assembly
  $$settings{order_seqs_by_name}=1;
  AKUtils::printSeqsToFile($newAssembly,$$settings{outfile},$settings);
  logmsg "Printed new assembly file $$settings{outfile}";

  return 1;
}

sub setupBuildDirectory{
  my($tempdir,$settings)=@_;
  system("rm -rf $$settings{tempdir}") if($$settings{force});
  warn "Warning: temporary directory already exists. -f to force a fresh start\n  $$settings{tempdir}" if(-e $$settings{tempdir});
  mkdir($$settings{tempdir});

  my $newSamPath="$tempdir/".basename($$settings{sam});
  my $newAssemblyPath="$tempdir/".basename($$settings{assembly});

  system("ln -s $$settings{sam} $newSamPath");
  system("ln -s $$settings{assembly} $newAssemblyPath");

  $$settings{sam}=$newSamPath;
  $$settings{assembly}=$newAssemblyPath;

  return 1;
}

# The database/reference genome file in the FASTA format was first indexed with the `index' 
sub indexAssembly{
  my($assembly,$settings)=@_;
  logmsg "Indexing assembly";
  if(!-e "$assembly.rsa"){
    my $command="bwa index -a bwtsw $assembly";
    command($command);
  } else {logmsg "Index exists; skipping";}

  return $assembly;
}

# the query illumina reads were aligned using the short read algorithm
sub queryIlluminaReads{
  my($settings)=@_;
  my $outfile="$$settings{tempdir}/aln_sa.sai";
  logmsg "Querying Illumina reads against assembly";
  if(-e $outfile){
    logmsg "Skipping because $outfile already exists. Use -f to force";
    return $outfile;
  }
  my $command="bwa aln $$settings{assembly} $$settings{illumina} > $outfile";
  command($command);
  
  return $outfile;
}

# The .sai output format was converted to the .sam output format as follows
sub convertSaiToSam{
  my($sai,$settings)=@_;
  my $outfile="$$settings{tempdir}/sam.aln";
  logmsg "Converting Sai to Sam format";
  if(-e $outfile){
    logmsg "Skipping because $outfile already exists. Use -f to force";
    return $outfile;
  }
  my $command="bwa samse $$settings{assembly} $sai $$settings{illumina} > $outfile";
  command($command);

  return $outfile;
}

# generate a .bam, .sorted.bam, .sorted.bam.bai and a database.fasta.fai file
sub createBam{
  my($sam,$settings)=@_;
  my $bamaln="$$settings{tempdir}/aln.bam";
  my $samindexed="$$settings{assembly}.fai";
  my $alnSortedPrefix="$$settings{tempdir}/aln.sorted";
  my $alnSorted="$alnSortedPrefix.bam"; # bam file
  my $bamIndex="$alnSorted.bai";

  logmsg "Converting SAM to BAM";

  # import the bam
  if(!-e $bamaln){
    command("samtools faidx $$settings{assembly}");
    command("samtools import $samindexed $sam $bamaln");
  } else {logmsg "Skipping bam import because file $bamaln already exists"; }

  # index the bam
  if(!-e $bamIndex){
    command("samtools sort $bamaln $alnSortedPrefix");
    command("samtools index $alnSorted");
  } else {logmsg "Skipping indexing the bam indexing because $alnSorted already exists";}

  return ($alnSorted,$bamIndex);
}

sub bamToFastq{
  my($bam,$bamIndex,$settings)=@_;

  my $out1="$$settings{tempdir}/mpileupout1.tmp";
  my $out2="$$settings{tempdir}/bcftoolsout2.tmp";
  my $fastqOutNonstandard="$$settings{tempdir}/outNonstandard.fastq";
  my $fastqOutStandard="$$settings{tempdir}/outStandard.fastq";
  my $fastqOut="$$settings{tempdir}/out.fastq";

  logmsg "Converting BAM to fastq";
  if(!-e $fastqOutNonstandard){
    #indexAssembly($$settings{assembly},$settings);
    # separate out these commands for debugging purposes
    command("samtools mpileup -uf $$settings{assembly} $bam > $out1");
    command("bcftools view -cg - < $out1 > $out2");
    command("vcfutils.pl vcf2fq < $out2 > $fastqOutNonstandard");
  } else {logmsg "$fastqOutNonstandard exists; skipping";}
  if(!-e $fastqOut || -s $fastqOut < 1){
    standardizeFastq($fastqOutNonstandard,$fastqOutStandard,$settings);
    splitFastqByGaps($fastqOutStandard,$fastqOut,$settings);
  }
  return $fastqOut;
}

# change a fastq with tons of newlines to the standard 4-line entries
sub standardizeFastq{
  my($fIn,$fOut,$settings)=@_;
  logmsg "Standardizing $fIn to $fOut";
  open(IN,$fIn) or die "Could not open $fIn because $!";
  open(OUT,">",$fOut) or die "Could not write to $fOut because $!";
  my $id;
  while(my $line=<IN>){
    if($line=~/^@/){
      my($entrySequence,$entryQuality)=("","");
      $id=$line;
      chomp($id);
      # grab the sequence
      my $linesOfSequence=0;
      while(<IN>){
        last if(/^\+\s*$/);
        $entrySequence.=$_;
        $linesOfSequence++;
      }
      # grab the quality
      my $linesOfQuality=0;
      while(<IN>){
        $entryQuality.=$_;
        last if(++$linesOfQuality == $linesOfSequence);
      }

      #cleanup of the entry
      $entrySequence=~s/\n//g;
      $entryQuality=~s/\n//g;
      print OUT "$id\n$entrySequence\n+\n$entryQuality\n";
    }
  }
  close OUT;
  close IN;
  return $fOut;
}

sub fastqToFastaQual{
  my($fastq,$settings)=@_;
  open(FASTQ,"<",$fastq) or die "Could not open $fastq because $!";
  my $i=0;
  my $seqId="";
  my $sequence="NNN";
  my($fasta,$qual);
  while(my $fastqEntry=<FASTQ>){
    $fastqEntry.=<FASTQ> for(1..3);
    my($fastaEntry,$qualEntry)=fastqEntryToFastaQual($fastqEntry,$settings);
    $qual.="$qualEntry\n";
    $fasta.="$fastaEntry\n";
  }
  close FASTQ;
  return($fasta,$qual);
}

sub fastqEntryToFastaQual{
  my($entry,$settings)=@_;
  $entry=~s/^\s+|\s+$//g;
  my($id,$sequence,undef,$qualIllumina)=split(/\s*\n\s*/,$entry);
  $id=~s/^@/>/;
  my @qual=split(//,$qualIllumina);
  @qual=map(ord($_)-33,@qual);
  for(my $i=60;$i<@qual;$i+=60){
    $qual[$i].="\n";
  }
  my $qualFasta=join(" ",@qual);
  
  $sequence=~s/(.{60})/$1\n/g;
  my $fastaEntry="$id\n$sequence";
  my $qualEntry="$id\n$qualFasta";
  return ($fastaEntry,$qualEntry) if wantarray;
  return $fastaEntry;
}

# split contigs by gaps
sub splitFastqByGaps{
  my($fastq,$fastqOut,$settings)=@_;
  $$settings{min_gap_size}||=6;
  my $contigBreak="N"x$$settings{min_gap_size};

  logmsg "Splitting contigs by gaps of size $$settings{min_gap_size} from $fastq to $fastqOut";
  open(IN,$fastq) or die "Could not open $fastq because $!";
  open(OUT,">",$fastqOut) or die "Could not write to $fastqOut because $!";
  while(my $seqId=<IN>){
    my $sequence=<IN>;
    <IN>; # + sign
    my $qual=<IN>;
    chomp($seqId,$sequence,$qual);
    
    # look for large enough gaps to split contigs
    my $seqLength=length($sequence);
    my (@gapStart,@gapStop);
    my $contigNum=0;
    for(my $i=0;$i<$seqLength;$i++){ 
      next if(uc(substr($sequence,$i,$$settings{min_gap_size})) ne $contigBreak);
      # step 1: the contig break has been found.
      $gapStart[$contigNum]=$i;
      # step 2: elongation
      $i+=$$settings{min_gap_size};
      for($i=$i;$i<$seqLength;$i++){
        last if(uc(substr($sequence,$i,1)) ne "N");
      }
      $gapStop[$contigNum]=$i-1;
      $contigNum++;
    }

    # extract contigs between gaps
    my @newContig;
    my $contigStart=0;
    for(my $i=0;$i<@gapStart;$i++){
      my $subseqId=$seqId."_subcontig$i";
      my $contigStop;
      if($i+1<@gapStart){
        $contigStop=$gapStart[$i+1]-1;
      } else { $contigStop=$seqLength; }
      if($contigStop<1){
        next;
      }
      $contigStart=$gapStop[$i]+1;
      my $newContig=substr($sequence,$contigStart,($contigStop-$contigStart+1));
      my $newQual=substr($qual,$contigStart,($contigStop-$contigStart+1));
      print OUT "$subseqId\n$newContig\n+\n$newQual\n";
    }

  }
  close OUT;
  close IN;
  return $fastqOut;
}

#########################
### utility methods
#########################

sub readPileupLine{
  my($line,$settings)=@_;
  my %x=(type=>'snp');
  my @line=split(/\t/,$line);
  ($x{'chr'},$x{'pos'},$x{refBase},$x{consensusBase},$x{consensusQuality},$x{snpQuality},$x{mappingQuality},$x{coverage},$x{bases},$x{quals})=@line;

  # indels
  if($x{refBase} eq '*'){
    $x{type}='indel';
    delete($x{bases});
    delete($x{quals});
    $x{consensusBase}=$line[8];
    $x{consensusBase}=$line[9] if($x{consensusBase} eq '*');
  }

  # if you want to look at each base
  if($$settings{pileup_splitBases}){
    $x{basesArr}=[split(//,$x{bases})];
    $x{qualsArr}=[split(//,$x{quals})];
    #TODO put in numeric values for quality values
  }
  
  return %x;
}
sub command{
  my ($command)=@_;
  logmsg "RUNNING COMMAND\n  $command";
  system($command);
  die "ERROR running command $command\n  With error code $?" if($?);
  return 1;
}

sub usage{
  "This script converts a sam to an assembly.
Usage: perl $0 -a reference.fasta -s assembly.sam -o assembly.fasta -q
  -s sam file
  -a assembly reference file file
  -o (optional) final output file
    default: $0.merged.fasta
  -t (optional) This is where temporary files will be stored. 
    Default: /tmp/....../

  No arguments should be given for the following options
  -f (optional) Force.
  -k (optional) keep temporary files around (about 2GB of files in my test run)
  -q to output quality files too. (assembly.fasta.qual)
  ";
}

local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };

