#!/usr/bin/env perl

# run_assembly_trimClean: trim and clean a set of raw reads
# Author: Lee Katz <lkatz@cdc.gov>
# TODO read .gz files
# TODO output .gz files

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

use Getopt::Long;
use File::Temp ('tempdir');
use File::Path;
use File::Spec;
use File::Copy;
use File::Basename;
use List::Util qw(min max sum shuffle);
#use Math::Round qw/nearest/;

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
  my $settings={
    numcpus=>getNumCPUs(),
    poly=>1, # number of reads per group (1=SE, 2=paired end)
    qualOffset=>33,
    # trimming
    min_quality=>35,
    bases_to_trim=>20, # max number of bases that will be trimmed on either side
    # cleaning
    min_avg_quality=>30,
    min_length=>62,# twice kmer length sounds good
  };
  
  GetOptions($settings,qw(poly=i infile=s outfile=s min_quality=i bases_to_trim=i min_avg_quality=i  min_length=i quieter notrim));
  
  my $infile=$$settings{infile} or die "Error: need an infile\n".usage($settings);
  my $outfile=$$settings{outfile} or die "Error: need an outfile\n".usage($settings);
  
  my @IOextensions=qw(.fastq .fastq.gz);
  my $outfileDir;
  ($$settings{outBasename},$outfileDir,$$settings{outSuffix}) = fileparse($outfile,@IOextensions);
  (undef,undef,$$settings{inSuffix}) = fileparse($infile,@IOextensions);
  my $singletonOutfile="$outfileDir/$$settings{outBasename}.singletons$$settings{outSuffix}";
  
  my $readsToCheck=5;
  my $warnings=checkFirstXReads($infile,$readsToCheck,$settings);
  warn "WARNING: There were $warnings warnings out of $readsToCheck" if($warnings);

  my $printQueue=Thread::Queue->new;
  my $singletonPrintQueue=Thread::Queue->new;
  my $printThread=threads->new(\&printWorker,$outfile,$printQueue,$settings);
  my $singletonPrintThread=threads->new(\&singletonPrintWorker,$singletonOutfile,$singletonPrintQueue,$settings);
  
  my $entryCount=qualityTrimFastqPoly($infile,$printQueue,$singletonPrintQueue,$settings);
  
  $singletonPrintQueue->enqueue(undef); # term signal
  $singletonPrintThread->join;
  $printQueue->enqueue(undef); # term signal
  $printThread->join;
  
  my $numGood=sum(values(%threadStatus));
  my $freq_isClean=$numGood/$entryCount;
  $freq_isClean=nearest(0.01,$freq_isClean);
  logmsg "Finished! $freq_isClean pass rate.";

  return 0;
}

# returns a reference to an array of fastq entries.
# These entries are each a string of the actual entry
sub qualityTrimFastqPoly($;$){
  my($fastq,$printQueue,$singletonQueue,$settings)=@_;

  logmsg "Trimming and cleaning a file $fastq with poly=$$settings{poly}";
  # initialize the threads
  my (@t,$t);
  my $Q=Thread::Queue->new;
  for(0..$$settings{numcpus}-1){
    $t[$t++]=threads->create(\&trimCleanPolyWorker,$Q,$printQueue,$singletonQueue,$settings);
  }

  # load all reads into the threads for analysis
  my $entryCount=0;
  my $linesPerGroup=4*$$settings{poly};
  my $moreLinesPerGroup=$linesPerGroup-1; # calculate that outside the loop to save CPU
  if($$settings{inSuffix}=~/\.fastq$/){
    open(FQ,'<',$fastq) or die "Could not open $fastq for reading: $!";
  }
  elsif($$settings{inSuffix}=~/\.fastq\.gz/){
    open(FQ,"gunzip -c $fastq | ") or die "Could not open $fastq for reading: $!";
  }
  else{
    die "Could not determine the file type for reading based on your extension $$settings{inSuffix}";
  }
  while(my $entry=<FQ>){
    $entry.=<FQ> for(1..$moreLinesPerGroup); # e.g. 8 lines total for paired end entry
    $Q->enqueue($entry);
    $entryCount++;
    if(!$$settings{quieter} && $entryCount%100000==0){
      my $numGood=sum(values(%threadStatus));
      my $freq_isClean=$numGood/$entryCount;
      $freq_isClean=nearest(0.01,$freq_isClean);
      logmsg "Finished loading $entryCount pairs or singletons ($freq_isClean pass rate)";
    }
  }
  close FQ;
  logmsg "Done loading: $entryCount entries loaded.";
  
  # stop the threads
  $Q->enqueue(undef) for(1..$$settings{numcpus});
  queue_status_updater($Q,$settings);

  # grab all the results from the threads
  for(@t){
    my $tmp=$_->join;
  }
  return $entryCount;
}

sub trimCleanPolyWorker{
  my($Q,$printQueue,$singletonQueue,$settings)=@_;
  my $tid="TID".threads->tid;
  my(@entryOut);
  ENTRY:
  while(defined(my $entry=$Q->dequeue)){
    my @read;
    my @entryLine=split(/\s*\n\s*/,$entry);
    for my $i (0..$$settings{poly}-1){
      ($read[$i]{id},$read[$i]{seq},undef,$read[$i]{qual})=splice(@entryLine,0,4);
    }
    
    # this is just my dumb way of making sure I trim by reference
    for my $read (@read){
      trimRead($read,$settings);
    }
    
    # see if these are signletons, and if so, add them to the singleton queue
    my @singleton;
    for(my $i=0;$i<$$settings{poly};$i++){
      # if one is not good, then just go through them all to find out which can be retained
      if(!read_is_good($read[$i],$settings)){
        for my $read (@read){
          next if(!read_is_good($read,$settings));
          my $singletonOut="$$read{id}\n$$read{seq}\n+\n$$read{qual}\n";
          $singletonQueue->enqueue($singletonOut);
          $threadStatus{$tid}+=1/$$settings{poly};
        }
        next ENTRY;
      }
    }
    
    $threadStatus{$tid}++;
    my $entryOut="";
    for my $read (@read){
      $entryOut.="$$read{id}\n$$read{seq}\n+\n$$read{qual}\n";
    }
    $printQueue->enqueue($entryOut);
  }
  #return \@entryOut;
}

sub printWorker{
  my($outfile,$Q,$settings)=@_;
  if($$settings{outSuffix}=~/\.fastq$/){
    open(OUT,">",$outfile) or die "Error: could not open $outfile for writing because $!";
  }
  elsif($$settings{outSuffix}=~/\.fastq\.gz/){
    open (OUT, "| gzip -c > $outfile") or die "Error: could not open $outfile for writing because $!";
  }
  else{
    die "Could not determine the file type for writing based on your extension $$settings{outSuffix}";
  }
  while(defined(my $line=$Q->dequeue)){
    print OUT $line;
  }
  close OUT;
  logmsg "Reads are in $outfile";
}
sub singletonPrintWorker{
  my($outfile,$Q,$settings)=@_;
  if($$settings{outSuffix}=~/\.fastq$/){
    open(OUT,">",$outfile) or die "Error: could not open $outfile for writing because $!";
  }
  elsif($$settings{outSuffix}=~/\.fastq\.gz/){
    open (OUT, "| gzip -c > $outfile") or die "Error: could not open $outfile for writing because $!";
  }
  else{
    die "Could not determine the file type for writing based on your extension $$settings{outSuffix}";
  }
  while(defined(my $line=$Q->dequeue)){
    print OUT $line;
  }
  close OUT;
  logmsg "Singleton reads are in $outfile";
}

# Trim a read using the trimming options
# The read is a hash of id, sequence, and qual
# Note: I think I mixed up 3 and 5?  Ugh.  Whatever.
sub trimRead{
  my($read,$settings)=@_;
  return if($$settings{notrim});
  my($numToTrim3,$numToTrim5,@qual)=(0,0,);
  @qual=map(ord($_)-$$settings{qualOffset},split(//,$$read{qual}));
  for(my $i=0;$i<$$settings{bases_to_trim};$i++){
    $numToTrim5++ if(sum(@qual[0..$i])/($i+1)<$$settings{min_quality});
    #last if($qual[$i]>=$$settings{min_quality});
  }
  #for(my $i=$$read{length}-$$settings{bases_to_trim};$i<@qual;$i++){
  for(my $i=0;$i<$$settings{bases_to_trim};$i++){
    my $pos=$$read{length}-$i;
    $numToTrim3++ if(sum(@qual[$pos..$#qual])/($i+1)<$$settings{min_quality});
    #last if($qual[$i]>=$$settings{min_quality});
  }
  if($numToTrim3 || $numToTrim5){
    #my $tmp= "TRIM  ==>$numToTrim5 ... $numToTrim3<==\n  >$$read{seq}<\n   $$read{qual}";
    $$read{seq}=substr($$read{seq},$numToTrim5,$$read{length}-$numToTrim5-$numToTrim3);
    $$read{qual}=substr($$read{qual},$numToTrim5,$$read{length}-$numToTrim5-$numToTrim3);
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
  my @qual=map(ord($_)-$$settings{qualOffset},split(//,$$read{qual}));
  my $avgQual=sum(@qual)/$readLength;
  return 0 if($avgQual<$$settings{min_avg_quality});
  return 1;
}

sub queue_status_updater{
  my($Q,$settings)=@_;
  STATUS_UPDATE:
  while(my $p=$Q->pending){
    logmsg "$p entries left to process";
    # shave off some time if the queue empties while we are sleeping
    for(1..5){
      sleep 2;
      last STATUS_UPDATE if(!$Q->pending);
    }
  }
  return 1;
}

sub checkFirstXReads{
  my($infile,$numReads,$settings)=@_;
  my $warning=0;
  # figure out if the min_length is larger than any of the first X reads.
  # If so, spit out a warning. Do not change the value. The user should be allowed to make a mistake.
  logmsg "Checking minimum sequence length param on the first 5 sequences.";
  if($$settings{inSuffix}=~/\.fastq$/){
    open(IN,'<',$infile) or die "Could not open $infile for reading: $!";
  }
  elsif($$settings{inSuffix}=~/\.fastq\.gz/){
    open(IN,"gunzip -c $infile | ") or die "Could not open $infile for reading: $!";
  }
  else{
    die "Could not determine the file type for reading based on your extension $$settings{inSuffix}";
  }
  my $i=1;
  while(my $line=<IN>){
    $line.=<IN> for (1..3);
    my ($id,$sequence)=(split("\n",$line))[0,1];
    my $length=length($sequence);
    if($length<$$settings{min_length}){
      warn "WARNING: Read $id has length greater than the minimum sequence length $$settings{min_length}. It is possible that all sequences will be filtered out.\n" if(!$$settings{quieter});
      $warning++;
    }
    last if($i++>=$numReads);
  }
  close IN;
  return $warning;
}

################
## utility
###############

sub nearest{
  return @_[0];
}

sub getNumCPUs() {
  my $num_cpus;
  open(IN, '<', '/proc/cpuinfo'); while (<IN>) { /processor\s*\:\s*\d+/ or next; $num_cpus++; } close IN;
  return $num_cpus || 1;
}

sub usage{
  my ($settings)=@_;
  "trim and clean a set of raw reads
  Usage: $0 -i reads.fastq -o reads.filteredCleaned.fastq [-p]
    -i input file in fastq format
    -o output file in fastq format
  Additional options
  
  -s to produce a singletons output file, for when one out of a paired end read is cleaned
  -p 1 or 2 (p for poly)
    1 for SE, 2 for paired end (PE) 
  -q for somewhat quiet mode (use 1>/dev/null for totally quiet)
  --notrim to skip trimming of the reads. Useful for assemblers that require equal read lengths.

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


