#! /usr/bin/env perl

die "Use nesoni instead";

# Adds Illumina reads to an assembly to make better base calls and improve the quality
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

  GetOptions($settings,qw(assembly=s illumina=s force tempdir=s keep outfile=s));

  # check for required parameters
  for my $param (qw(assembly illumina)){
    $$settings{$param} || die "Error: need $param parameter\n".usage();
  }
  # make some filenames absolute
  for my $param (qw(assembly illumina)){
    $$settings{$param} = File::Spec->rel2abs($$settings{$param});
  }

  # set up the directory structure
  $$settings{tempdir}  ||= tempdir(File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($$settings{keep}));
  $$settings{tempdir}=File::Spec->rel2abs($$settings{tempdir});
  $$settings{outfile}||="$0.merged.fasta";
  logmsg "Temporary directory is $$settings{tempdir} and final output will be in $$settings{ChangeDir}";

  # start the computation

  setupBuildDirectory($$settings{tempdir},$settings);
  indexAssembly($settings);

  my $sai=queryIlluminaReads($settings);

  my $sam=convertSaiToSam($sai,$settings);

  my ($bam,$bamIndex)=createBam($sam,$settings);

  my $snpcalls=predictSnps($bam,$bamIndex,$settings);

  my %baseChanges=makeNewBasecalls($snpcalls,$settings);
  my $newAssembly=putChangesIntoAssembly(\%baseChanges,$settings);

  #TODO print out the new assembly
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

  my $newAssPath="$tempdir/".basename($$settings{assembly});
  my $newIllPath="$tempdir/".basename($$settings{illumina});

  system("ln -s $$settings{assembly} $newAssPath");
  system("ln -s $$settings{illumina} $newIllPath");

  $$settings{assembly}=$newAssPath;
  $$settings{illumina}=$newIllPath;

  return 1;
}

# The database/reference genome file in the FASTA format was first indexed with the `index' 
sub indexAssembly{
  my($settings)=@_;
  logmsg "Indexing assembly";
  if(!-e "$$settings{assembly}.rsa"){
    my $command="bwa index -a bwtsw $$settings{assembly}";
    command($command);
  } else {logmsg "Index exists; skipping";}

  return $$settings{assembly};
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

# generate the pileup output of SNP/indel predictions
sub predictSnps{
  my($bam,$bamIndex,$settings)=@_;
  my $rawSnps="$$settings{tempdir}/sam.raw.snps";
  my $filteredSnps="$$settings{tempdir}/sam.filtered.snps";
  my $pileup="$$settings{tempdir}/sam.pileup.snps";

  logmsg "Creating pileup file";
  #print "BAM FILE IS $bam and ASSEMBLY FILE IS $$settings{assembly}\n";die;
  # samtools mpileup -ugf ref.fa aln1.bam aln2.bam | bcftools view -bvcg - > var.raw.bcf
  # bcftools view var.raw.bcf | vcfutils.pl varFilter -D100 > var.flt.vcf
  # works -- samtools mpileup -ugf /netapp_beryl_users/gzu2/projects/assembly/tmp/2010V-1031.assembly.fasta /netapp_beryl_users/gzu2/projects/assembly/tmp/2010V-1031.assembly.fasta.bam.bai | bcftools view -bvcg - > var.raw.bcf
  if(!-e $pileup){
    command("samtools pileup -vcf $$settings{assembly} $bam > $rawSnps");
    command("samtools.pl varFilter -D100 < $rawSnps > $filteredSnps") if(!-e $pileup);
    command("awk '(\$3==\"*\"&&\$6>=50) || (\$3!=\"*\"&&\$6>=20)' $filteredSnps > $pileup");
  } else {logmsg "$pileup exists; skipping";}

  return $pileup;
}

# determine the needed coverage before a base can even be considered to be changed
sub coverageThreshold{
  my ($snpcalls,$settings)=@_;

  # just start off on 20 for now, using Eishita's research
  return 20;
}

sub percentNeededForSnpCall{
  my($snpcalls,$minCoverage,$settings)=@_;
  
  # just start off on 90 for now, using Eishita's research
  return 90;
}

# list where all base call changes are
sub makeNewBasecalls{
  my($snpcalls,$settings)=@_;
  my $minCoverage=coverageThreshold($snpcalls,$settings);
  my $percentNeededForSnpCall=percentNeededForSnpCall($snpcalls,$minCoverage,$settings);

  my %changes=();

  my $seqs=AKUtils::readMfa($$settings{assembly});
  open(SNPS,$snpcalls) or die "Cannot open pileup snpcalls file $snpcalls beacuse $!";
  while(my $line=<SNPS>){
    my %line=readPileupLine($line,$settings);
    #print Dumper \%line;die;
    next if($line{coverage}<$minCoverage);
    
    # mark bases for change, if they meet the requirements
    #TODO remove or ignore bases with low quality values

    # how many reads need to support you for a base change?
    my $supportNeededForChange=($percentNeededForSnpCall/100)*$line{coverage};

    my $consensusBase=$line{consensusBase};
    #my $poskey=join("_",$line{'chr'},$line{'pos'}); # position is defined by contig/position
    my($chr,$pos)=($line{'chr'},$line{'pos'});
    if($line{type} eq 'snp'){
      my $numBasesThatSupportTheSnp=($line{bases}=~s/$consensusBase//gi);
      if($numBasesThatSupportTheSnp>$supportNeededForChange){
        $changes{$chr}{$pos}=$line{consensusBase};
      }
    } elsif($line{type} eq 'indel'){
      $changes{$chr}{$pos}=$line{consensusBase};
    } else { die "Unknown type '$line{type}' in the pileup file"};
  }
  close SNPS;

  return %changes;
}

# make a new assembly based on these changes
sub putChangesIntoAssembly{
  my($baseChanges,$settings)=@_;
  my $assembly=AKUtils::readMfa($$settings{assembly});
  my $totalSubs=0;

  # make each change
  # changes are 1-based and the assembly is 0-based (b/c it's a regular string)
  while(my($contig,$changehash)=each(%$baseChanges)){
    my($contigSubs)=(0);
    my $sequence=[split(//,$$assembly{$contig})]; # make changes in an array

    # make changes in reverse numerical order so that indels don't mess up the coordinates
    # of any subsequent changes
    my @sortedPos=sort{
      #local($a,$b); # don't alter the original variables
      my($A,$B)=($a,$b);
      $A=~s/\+|\-//;
      $B=~s/\+|\-//;
      $B<=>$A
    } keys(%$changehash);
    for my $pos(@sortedPos){
      my $change=$$changehash{$pos};
      if($change!~/\+|\-/){
        makeSubstitution($sequence,$pos-1,$change,$settings);
      } elsif($change=~/\+/){
        makeInsertion($sequence,$pos-1,$change,$settings);
      } elsif($change=~/\-/){
        makeDeletion($sequence,$pos-1,$change,$settings);
      } 
      $contigSubs++;
    }
    $$assembly{$contig}=join("",@$sequence);

    logmsg "Made $contigSubs substitutions and indels for contig $contig";
    $totalSubs+=$contigSubs;
  }
  logmsg "Finished making substitutions and indels. $totalSubs total";

  return $assembly;
}


#########################
### utility methods
#########################

sub makeDeletion{
  my($sequence,$position,$nts,$settings)=@_;
  $nts=~s/\+|\-//;
  $position++; # 1-based
  # check to make sure the nt is what it says it should be
  if(substr($nts,0,1) ne $$sequence[$position]){
    print "There is a mismatch in the deletion at pos $position where it is $$sequence[$position] but the pileup file says ".substr($nts,0,1)."\n";
  }
  splice(@$sequence,$position,length($nts)); # remove nts starting at the ref position (?)
}
sub makeInsertion{
  my($sequence,$position,$nts,$settings)=@_;
  $nts=~s/\+|\-//;
  splice(@$sequence,$position+1,0,$nts); # put the nts into the position after the ref position
}

# make a change in the sequence
# 0-based coordinates
sub makeSubstitution{
  my($sequence,$position,$nt,$settings)=@_;
  $$sequence[$position]=$nt;
}

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
  "This script adds illumina reads to an existing assembly.
Usage: perl $0 -a assembly.fasta -i illuminaReads
  -a assembly to improve, in fasta format
  -i illumina reads, in fastq format
  -o (optional) final output file
    default: $0.merged.fasta
  -t (optional) This is where temporary files will be stored. 
    Default: /tmp/$0/

  No arguments should be given for the following options
  -f (optional) Force.
  -k (optional) keep temporary files around (about 2GB of files)
  ";
}

local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };

