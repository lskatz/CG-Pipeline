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

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {print STDOUT "$0: ".(caller(1))[3].": @_\n";}

exit(main());

sub main() {
  my $settings={};
  GetOptions($settings,qw(pairedEnd infile=s outfile=s min_quality=i bases_to_trim=i min_avg_quality=i min_avg_quality_window=i min_length=i));
  $$settings{numcpus}||=AKUtils::getNumCPUs();

  # trimming
  $$settings{min_quality}||=35;
  $$settings{bases_to_trim}||=10; # max number of bases that will be trimmed on either side
  # cleaning
  $$settings{min_avg_quality}||=30;
  $$settings{min_avg_quality_window}||=70;
  $$settings{min_length}||=62; # twice kmer length sounds good
  
  
  my $infile=$$settings{infile} or die "Error: need an infile\n".usage();
  my $outfile=$$settings{outfile} or die "Error: need an outfile\n".usage();

  logmsg "Trimming and cleaning a read file $infile";
  my $seqs=qualityTrimFastq($infile,$settings);
  
  open(OUT,">",$outfile);
  print OUT $_ for(@$seqs);
  close OUT;
  
  logmsg "Reads are in $outfile";

  return 0;
}

# returns a reference to an array of fastq entries.
# These entries are each a string of the actual entry
sub qualityTrimFastq($;$){
  my($fastq,$settings)=@_;

  # initialize the threads
  my (@t,$t);
  my $Q=Thread::Queue->new;
  for(0..$$settings{numcpus}-1){
    if($$settings{pairedEnd}){
      $t[$t++]=threads->create(\&trimCleanPEWorker,$Q,$settings);
    } else {
      die "Error: non paired end does not work right now. Use -p";
    }
  }

  # load all reads into the threads for analysis
  my $entryCount=0;
  open(FQ,'<',$fastq) or die "Could not open $fastq for reading: $!";
  while(my $entry=<FQ>){
    $entry.=<FQ> for(1..7); # 8 lines total for paired end entry
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

  my $numCleaned=$entryCount-@entry;
  logmsg "$numCleaned entries out of $entryCount removed due to low quality";
  return \@entry;
}

sub trimCleanPEWorker{
  my($Q,$settings)=@_;
  logmsg "Launching thread TID".threads->tid;
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
    push(@entryOut,$entry);
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
  $$read{seq}=substr($$read{seq},$numToTrim3,$$read{length}-$numToTrim3-$numToTrim5);
  $$read{qual}=substr($$read{qual},$numToTrim3,$$read{length}-$numToTrim3-$numToTrim5);

  return $read;
}



# Deterimine if a single read is good or bad (0 or 1)
# judging by the cleaning options
sub read_is_good{
  my($read,$settings)=@_;
  return 0 if(length($$read{seq})<$$settings{min_length});
  
  # avg quality
  my $qual_is_good=1;
  my @qual=map(ord($_)-33,split(//,$$read{qual}));
  my $length=scalar(@qual)-$$settings{min_avg_quality_window};
  for(my $i=0;$i<$length;$i++){
    my @subQual=@qual[$i..($$settings{min_avg_quality_window}+$i-1)];
    my $avgQual=sum(@subQual)/$$settings{min_avg_quality_window};
    if($avgQual<$$settings{min_avg_quality}){
      #print "$i: ".join(".",@subQual)."\n$avgQual <=> $$settings{min_avg_quality}\n";die;
      return 0;
    }
  }
  return 1;
}


sub usage{
  my ($settings)=@_;
  "trim and clean a set of raw reads
  Usage: $0 -i reads.fastq -o reads.filteredCleaned.fastq [-p]
    -i input file in fastq format
    -o output file in fastq format
  Additional options
  
  -p Use this switch if the reads are paired-end (currently on by default)

  Use phred scores (e.g. 20 or 30) or length in base pairs if it says P or L, respectively
  --min_quality P             # trimming
  --bases_to_trim L           # trimming
  --min_avg_quality P         # cleaning
  --min_avg_quality_window L  # cleaning
  --min_length L              # cleaning
  "
}

