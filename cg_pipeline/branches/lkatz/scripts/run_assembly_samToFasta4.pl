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

sub bamToFastq{
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
    next if(!@$bases);
    my $consensus=join("",@$bases);
    $consensus=~s/(.{60})/$1\n/g;
    print FASTAOUT ">$rseq\n$consensus\n";
  }
  close FASTAOUT;
  logmsg "Fasta is now in $fastaOut";die;

  splitFastqByGaps($fastqOutStandard,$fastqOut,$settings);
  return $fastqOut;
}

sub mpileupWorker{
  my($Q,$mpileupPrintQueue,$minDepth,$maxDepth,$settings)=@_;
  my $tid="TID".threads->tid;
  my $tidNum=threads->tid;
  while(defined(my $samLine=$Q->dequeue)){
    chomp($samLine);
    my (%baseCount);
    my($rseq,$pos,undef,$depth,$pileup,$quality)=split(/\t/,$samLine);
    $depth=1 if($depth<1); # avoid divide by 0 error

    #die "TODO consider insertions";
    my ($p,$q)=_pileupStrToArr($pileup,$quality,$settings);
    $baseCount{$_}++ for (@$p);
    $baseCount{$_}/=$depth for(keys(%baseCount));
    my $topPossibility=(sort({$baseCount{$b}<=>$baseCount{$a}} keys(%baseCount)))[0];
    my $topPossibilityPercent=$baseCount{$topPossibility};

    my $baseCall="";
    my $qualityCall="";
    if($depth<$minDepth || $depth > $maxDepth){
      $baseCall="N";
      $qualityCall=0;
    } elsif($topPossibility=~/^([+-])/){  # TODO take care of insertions
      die "This is an indel site\n$pileup\n";
      if($1 eq '+'){
      } else { # eq '-'
        
      }
    } else{
      if($topPossibilityPercent>$$settings{pileup_min_frequency}){
        $baseCall=$topPossibility;
        $qualityCall=_qualityCall($baseCall,$p,$q,$settings);
        if($baseCall=~/[+-]/){
          die "!!! Need to understand how to deal with insertions. Found this: $baseCall";
        }
      } else {
        $baseCall="N";
        $qualityCall=0;
      }
    }

    $mpileupPrintQueue->enqueue(join("\t",$rseq,$pos,$baseCall,$qualityCall)."\n");
  }
  logmsg "$tid DONE";
  return 1;
}

# calculate the quality of an NT, based off of the pileup
sub _qualityCall{
  my($nt,$p,$q,$settings)=@_;
  my $qual=1; # start with 100% and go from there

  # calculate the quality score of all correct calls combined
  # formula: 1 minus product of all error percentages
  # subtract off the mismatch quality scores
  my $numBases=@$p;
  for(my $i=0;$i<$numBases;$i++){
    # determine quality
    if($$p[$i] eq $nt){
      $qual*=$$q[$i];
    } else {
      $qual/=$$q[$i];
    }
  }
  $qual=int(-10*log($qual||0.0000000001));
  #print Dumper [[join(" ",@p)],[join(" ",@q)],$qual];die;
  $qual=60 if($qual>60);
  #convert back to letter code
  $qual=chr($qual+33);
  return $qual;
}

# Converts a pileup string and quality string to arrays.
# Ignores forward/reverse distinctions.
sub _pileupStrToArr{
  my($pileup,$quality,$settings)=@_;
  my $qual=1; # start with 100% and go from there
  $pileup=uc($pileup);
  my @p=split(//,$pileup);
  my @q=split(//,$quality);
  @q=map(10**(-1*(ord($_)-33)/10),@q);
  my $numBases=@p;

  my (@outQ,@outP);
  my $j=0;
  for(my $i=0;$i<$numBases;$i++){
    my $p=$p[$i];
    my $q=$q[$j] || 1;
    #die "Internal error\n". Dumper("i $i ($p)","j $j ($q)","len p ".length($pileup),"len q ".length($quality),$pileup,$quality) if(!$q[$j] || $q[$j]<=0);

    ## take care of special characters
    # indels
    if($p=~/([+-])/){
      my $insLength=$p[$i+1]; # TODO WARNING assumes single-digit insertion
      $p.=join("",@p[($i+1)..($i+1+$insLength)]); # $p is already + or -
      $i+=length($p)-2; # $insLength+1; # increment the number of bases inserted plus the number
    }
    # start/stops
    elsif($p=~/([\^\$])/){
      my $mappingQual="";
      $mappingQual=$p[++$i] if($1 eq '^');
      $p.=$mappingQual;
    }
    next if(!defined($p) || $p=~/^\s*$/);
    push(@outP,$p);
    push(@outQ,$q);

    $j++; # increment index for qual scores
  }

  die "INTERNAL ERROR: pileup array length does not match quality scores array length\n".join(" ",@outP[0..12])."\n".join(" ",@outQ[0..12])."\n" if(scalar(@outP) != scalar(@outQ));

  return (\@outP,\@outQ);
}

sub fastqToFastaQual{
  my($fastq,$settings)=@_;

  open(FASTQ,"<",$fastq) or die "Could not open $fastq because $!";
  my $i=0;
  my $seqId="";
  my($fasta,$qual);
  while(my $fastqEntry=<FASTQ>){
    $fastqEntry.=<FASTQ> for(1..3);
    my($fastaEntry,$qualEntry)=fastqEntryToFastaQual($fastqEntry,$settings);
    next if(!$fastaEntry);
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
  return 0 if(!$qualIllumina || length($sequence) < 50); # return false if there's nothing here
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
  die "ERROR: need min_gap_size" if(!$$settings{min_gap_size});

  my $contigBreak="N"x$$settings{min_gap_size};
  my $qualContigBreak=chr(33)x$$settings{min_gap_size};

  logmsg "Splitting contigs by gaps of size $$settings{min_gap_size} from $fastq to $fastqOut";
  open(IN,$fastq) or die "Could not open $fastq because $!";
  open(OUT,">",$fastqOut) or die "Could not write to $fastqOut because $!";
  while(my $seqId=<IN>){
    my $sequence=<IN>;
    <IN>; # + sign
    my $qual=<IN>;
    chomp($seqId,$sequence,$qual);

    # Lowercase nts are low quality in either depth or mapping quality, according to samtools. 
    # Turn them into Ns if the user requests it.
    $sequence=~tr/atcgn/N/ if($$settings{mask});
    
    # look for large enough gaps to split contigs
    # A hack: start the sequence with the gap string
    $sequence="$contigBreak$sequence";
    $qual="$qualContigBreak$qual";
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

