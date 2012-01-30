#!/usr/bin/env perl

# run_assembly_trimClean: trim and clean a set of raw reads
# Author: Lee Katz <lkatz@cdc.gov>

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
use AKUtils;

use Getopt::Long;
use File::Temp ('tempdir');
use File::Path;
use File::Spec;
use File::Copy;
use File::Basename;
use List::Util qw(min max sum shuffle);
use CGPipelineUtils;
use Data::Dumper;
use threads;
use Thread::Queue;
use threads::shared;

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {print STDOUT "$0: ".(caller(1))[3].": @_\n";}

my %threadStatus:shared;

exit(main());

sub main() {
  my $settings={};
  GetOptions($settings,qw(pairedEnd infile=s outfile=s min_quality=i bases_to_trim=i min_avg_quality=i  min_length=i));
  $$settings{numcpus}||=AKUtils::getNumCPUs();

  # trimming
  $$settings{min_quality}||=35;
  $$settings{bases_to_trim}||=10; # max number of bases that will be trimmed on either side
  # cleaning
  $$settings{min_avg_quality}||=30;
  $$settings{min_length}||=62; # twice kmer length sounds good
  
  
  my $infile=$$settings{infile} or die "Error: need an infile\n".usage($settings);
  my $outfile=$$settings{outfile} or die "Error: need an outfile\n".usage($settings);

  my $seqs=[];
  $seqs=qualityTrimFastqPE($infile,$settings) if($$settings{pairedEnd});
  $seqs=qualityTrimFastqSE($infile,$settings) if(!$$settings{pairedEnd});
  
  open(OUT,">",$outfile) or die "Error: Could not open $outfile for writing";
  print OUT $_ for(@$seqs);
  close OUT;
  logmsg "Reads are in $outfile";

  return 0;
}

# returns a reference to an array of fastq entries.
# These entries are each a string of the actual entry
sub qualityTrimFastqPE($;$){
  my($fastq,$settings)=@_;

  logmsg "Trimming and cleaning a paired end read file $fastq";
  # initialize the threads
  my (@t,$t);
  my $Q=Thread::Queue->new;
  for(0..$$settings{numcpus}-1){
    $t[$t++]=threads->create(\&trimCleanPEWorker,$Q,$settings);
  }

  # load all reads into the threads for analysis
  my $entryCount=0;
  open(FQ,'<',$fastq) or die "Could not open $fastq for reading: $!";
  while(my $entry=<FQ>){
    $entry.=<FQ> for(1..7); # 8 lines total for paired end entry
    $Q->enqueue($entry);
    $entryCount++;
    if($entryCount%100000==0){
      my $numGood=sum(values(%threadStatus));
      my $freq_isClean=$numGood/$entryCount;
      logmsg "Finished loading $entryCount pairs ($freq_isClean pass rate)";
    }
  }
  close FQ;
  logmsg "Done loading: $entryCount entries loaded.";
  
  # stop the threads
  $Q->enqueue(undef) for(1..$$settings{numcpus});
  # status update
  STATUS_UPDATE:
  while(my $p=$Q->pending){
    logmsg "$p entries left to process";
    # shave off some time if the queue empties while we are sleeping
    for(1..5){
      sleep 2;
      last STATUS_UPDATE if(!$Q->pending);
    }
  }

  # grab all the results from the threads
  my @entry;
  for(@t){
    my $tmp=$_->join;
    push(@entry,@$tmp);
  }

  my $numEntriesLeft=@entry;
  my $numCleaned=$entryCount-$numEntriesLeft;
  logmsg "$numCleaned entries out of $entryCount removed due to low quality ($numEntriesLeft left)";
  return \@entry;
}
sub qualityTrimFastqSE($;$){
  my($fastq,$settings)=@_;

  logmsg "Trimming and cleaning a single end read file $fastq";
  # initialize the threads
  my (@t,$t);
  my $Q=Thread::Queue->new;
  for(0..$$settings{numcpus}-1){
    $t[$t++]=threads->create(\&trimCleanSEWorker,$Q,$settings);
  }
  # load all reads into the threads for analysis
  my $entryCount=0;
  open(FQ,'<',$fastq) or die "Could not open $fastq for reading: $!";
  while(my $entry=<FQ>){
    $entry.=<FQ> for(1..3); # 4 lines total for single end entry
    $Q->enqueue($entry);
    $entryCount++;
    logmsg "Finished loading $entryCount pairs" if($entryCount%100000==0);
  }
  close FQ;
  logmsg "Done loading: $entryCount entries loaded.";
  
  # stop the threads
  $Q->enqueue(undef) for(1..$$settings{numcpus});
  # status update
  STATUS_UPDATE:
  while(my $p=$Q->pending){
    logmsg "$p entries left to process";
    # shave off some time if the queue empties while we are sleeping
    for(1..5){
      sleep 2;
      last STATUS_UPDATE if(!$Q->pending);
    }
  }

  # grab all the results from the threads
  my @entry;
  for(@t){
    my $tmp=$_->join;
    push(@entry,@$tmp);
  }

  my $numEntriesLeft=@entry;
  my $numCleaned=$entryCount-$numEntriesLeft;
  logmsg "$numCleaned entries out of $entryCount removed due to low quality ($numEntriesLeft left)";
  return \@entry;
}

sub trimCleanPEWorker{
  my($Q,$settings)=@_;
  my $tid="TID".threads->tid;
  logmsg "Launching thread $tid";
  my(@entryOut);
  while(defined(my $entry=$Q->dequeue)){
    my($id1,$seq1,undef,$qual1,$id2,$seq2,undef,$qual2)=split(/\s*\n\s*/,$entry);
    my %read1=(id=>$id1,seq=>$seq1,qual=>$qual1,length=>length($seq1));
    my %read2=(id=>$id2,seq=>$seq2,qual=>$qual2,length=>length($seq2));
    trimRead(\%read1,$settings);
    trimRead(\%read2,$settings);

    #cleaning stage
    # TODO allow singletons to pass
    next if(!read_is_good(\%read1,$settings));
    next if(!read_is_good(\%read2,$settings));

    my $entryOut="$read1{id}\n$read1{seq}\n+\n$read1{qual}\n$read2{id}\n$read2{seq}\n+\n$read2{qual}\n";
    $threadStatus{$tid}++;
    push(@entryOut,$entryOut);
  }
  return \@entryOut;
}
sub trimCleanSEWorker{
  my($Q,$settings)=@_;
  logmsg "Launching thread TID".threads->tid;
  my(@entryOut);
  while(defined(my $entry=$Q->dequeue)){
    my($id1,$seq1,undef,$qual1)=split(/\s*\n\s*/,$entry);
    my %read1=(id=>$id1,seq=>$seq1,qual=>$qual1,length=>length($seq1));
    trimRead(\%read1,$settings);

    #cleaning stage
    next if(!read_is_good(\%read1,$settings));

    my $entryOut="$read1{id}\n$read1{seq}\n+\n$read1{qual}\n";
    push(@entryOut,$entryOut);
  }
  return \@entryOut;
}

# Trim a read using the trimming options
# The read is a hash of id, sequence, and qual
sub trimRead{
  my($read,$settings)=@_;
  my($numToTrim3,$numToTrim5,@qual);
  $numToTrim3=$numToTrim5=0;
  @qual=map(ord($_)-33,split(//,$$read{qual}));
  for(my $i=0;$i<$$settings{bases_to_trim};$i++){
    $numToTrim3++ if($qual[$i]<$$settings{min_quality});
    last if($qual[$i]>=$$settings{min_quality});
  }
  for(my $i=$$read{length}-$$settings{bases_to_trim};$i<@qual;$i++){
    $numToTrim5++ if($qual[$i]<$$settings{min_quality});
    last if($qual[$i]>=$$settings{min_quality});
  }
  if($numToTrim3 || $numToTrim5){
    #my $tmp= "TRIM  ==>$numToTrim3 ... <==$numToTrim5\n  >$$read{seq}<";
    $$read{seq}=substr($$read{seq},$numToTrim3,$$read{length}-$numToTrim3-$numToTrim5);
    $$read{qual}=substr($$read{qual},$numToTrim3,$$read{length}-$numToTrim3-$numToTrim5);
    #print "$tmp\n  >$$read{seq}<\n\n";
  }

  return %$read;
}



# Deterimine if a single read is good or bad (0 or 1)
# judging by the cleaning options
sub read_is_good{
  my($read,$settings)=@_;
  my $readLength=length($$read{seq});
  #die "short read: $$read{seq}\n$readLength < $$settings{min_length}\n" if($readLength<$$settings{min_length});
  return 0 if($readLength<$$settings{min_length});
  
  # avg quality
  my $qual_is_good=1;
  my @qual=map(ord($_)-33,split(//,$$read{qual}));
  my $avgQual=sum(@qual)/$readLength;
  return 0 if($avgQual<$$settings{min_avg_quality});
  return 1;
}


sub usage{
  my ($settings)=@_;
  "trim and clean a set of raw reads
  Usage: $0 -i reads.fastq -o reads.filteredCleaned.fastq [-p]
    -i input file in fastq format
    -o output file in fastq format
  Additional options
  
  -p Use this switch if the reads are paired-end 

  Use phred scores (e.g. 20 or 30) or length in base pairs if it says P or L, respectively
  --min_quality P             # trimming
    default: $$settings{min_quality}
  --bases_to_trim L           # trimming
    default: $$settings{bases_to_trim}
  --min_avg_quality P         # cleaning
    default: $$settings{min_avg_quality}
  --min_length L              # cleaning
    default: $$settings{min_length}
  "
}

