#!/usr/bin/env perl

# Predict CRISPRs in a fasta file
# Author: Lee Katz <lkatz@cdc.gov>

# TODO multithread
# TODO validate scoring mechanism
# TODO make GFF source tag for each contig

use strict;
use warnings;
use FindBin;
use lib "$FindBin::RealBin/../lib";
$ENV{PATH} = "$FindBin::RealBin:".$ENV{PATH};
use AKUtils qw/logmsg/;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use File::Copy qw/copy/;
use Bio::AlignIO;
use POSIX qw/mkfifo/;

use threads;
use Thread::Queue;

local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
exit(main());

sub main{
  local $0=fileparse $0;
  my $settings={
    appname=>"cgpipeline",
    minCrisprSize=>23,
    maxCrisprSize=>55,
  };
  GetOptions($settings,qw(outfile=s help windowSize=i stepSize=i tempdir=s minCrisprSize=i maxCrisprSize=i));
  $$settings{windowSize}||=$$settings{maxCrisprSize}*3;
  $$settings{stepSize}||=$$settings{windowSize}-2*$$settings{maxCrisprSize};
  $$settings{numcpus}||=AKUtils::getNumCPUs();
  $$settings{repeatMatch}=AKUtils::fullPathToExec("repeat-match");
  die usage($settings) if(@ARGV<1 || $$settings{help});
  my $outfile=$$settings{outfile} || "$0.gff";
  my $infile=$ARGV[0];
  die "ERROR: Could not find infile $infile" if(!-e $infile);
  $$settings{tempdir}||=AKUtils::mktempdir();
  logmsg "Temporary directory is $$settings{tempdir}";
  logmsg "Detecting CRISPRs in $infile";

  my $contigs=AKUtils::readMfa($infile,{first_word_only=>1});

  # this is the slow step. Running a repeat finder.
  my $smallTandemRepeatSeqs=findSmallTandemRepeats($contigs,$settings);

  # Finding the actual identified repeats is fast
  my $positions=findSequencePositionsInGenome($smallTandemRepeatSeqs,$contigs,$settings);
  $positions=removeOverlappingPositions($positions,$settings);
  my $gff=positionsToGff($positions,$contigs,$settings);
  copy($gff,$outfile);
  logmsg "Finished. GFF is in $outfile";
  return 0;
}

sub removeOverlappingPositions{
  my($positions,$settings)=@_;
  my %nonOverlapping;
  my %crisprStartStop;
  while(my($sequence,$drArr)=each(%$positions)){
    my @dr=sort({return 0 if($$a[0] ne $$b[0]); $$a[1] <=> $$b[1]} @$drArr);
    my $numDr=@dr;
    my $lo=$dr[0][1];
    my $hi=$dr[$numDr-1][2];

    # if this lo/hi is in another crispr, then don't use it
    my $isGood=1;
    for(values(%crisprStartStop)){
      my($contig,$start,$stop)=@$_;
      # inside
      $isGood=0 if($lo>=$start && $hi<=$stop);
      $isGood=0 if($lo<=$stop && $lo=>$start);
      $isGood=0 if($hi<=$stop && $hi=>$start);
    }
    next if(!$isGood);
    
    # PASSED: include this in the positions hash
    $crisprStartStop{$sequence}=[$dr[0][0],$lo,$hi];
    $nonOverlapping{$sequence}=$drArr;
  }

  return \%nonOverlapping;
}

sub filterRepeatSeqs{
  my($seqs,$settings)=@_;
  my %badSeqIndex;
  my $numSeqs=@$seqs;
  for(my $i=0;$i<$numSeqs;$i++){
    my $seqI=$$seqs[$i];
    for(my $j=$i+1;$j<$numSeqs;$j++){
      my $seqJ=$$seqs[$j];
      if($seqI=~/$seqJ/){ # if $seqJ is inside of $seqI, then $seqJ is a more likely, smaller repeat
        $badSeqIndex{$i}=1;
      }
    }
  }
  
  my @filteredSeq;
  for(my $i=0;$i<$numSeqs;$i++){
    next if($badSeqIndex{$i});
    push(@filteredSeq,$$seqs[$i]);
  }
  return \@filteredSeq;
}

sub positionsToGff{
  my($positions,$contigs,$settings)=@_;
  my $gff="$$settings{tempdir}/crispr.gff";
  open(GFF,">",$gff) or die "Could not open temporary GFF:$!";
  print GFF "##gff-version   3\n";
  print GFF "##feature-ontology http://song.cvs.sourceforge.net/viewvc/*checkout*/song/ontology/so.obo?revision=1.263\n";
  while(my($id,$sequence)=each(%$contigs)){
    my $start=1;
    my $end=length($sequence);
    print GFF "##sequence-region $id $start $end\n";
  }
  print GFF "###\n"; # indicates the end of directives
  ## done with GFF directives
  ###########################

  # contigs
  while(my($id,$sequence)=each(%$contigs)){
    my $start=1;
    my $end=length($sequence);
    print GFF join("\t",$id,".","contig",$start,$end,".","+",".","ID=$id")."\n";
  }
  # CRISPRs
  my $crisprCounter=1;
  my $drCounter=1;
  my $spacerCounter=1;
  for my $drSeq(keys(%$positions)){
    my $crispr=$$positions{$drSeq};
    $crispr=[sort({$$a[1] <=> $$b[1]} @$crispr)];
    my $numDrs=@$crispr;
    next if($numDrs<2); # need 2+ DRs for a CRISPR

    # formulate the overall crispr element score. This needs to be tested.
    my $score=1;
    my $crisprId=sprintf("crispr%s",$crisprCounter);
    my ($crisprStart,$crisprStop)=($$crispr[0][1]+1,$$crispr[-1][2]+1);
    $score=$score*.98 if($numDrs<5); # less confidence for few repeats
    $score=$score*.95 if($numDrs<4); # less confidence for few repeats
    $score=$score*.95 if($numDrs<3); # less confidence for few repeats
    my $alnScore=drAlignmentScore($crispr,$contigs,$settings);
    $score=$score*$alnScore;
    my $drLength=$$crispr[0][2]-$$crispr[0][1];
    # penalize for each out-of whack DR/spacer ratio
    for(my $i=0;$i<@$crispr-1;$i++){
      my $spacerLength=$$crispr[$i+1][1]-1-$$crispr[$i][2]+1;
      $score=$score*.6 if($spacerLength>2.5*$drLength || $spacerLength<0.6*$drLength);
    }

    # don't accept low-quality CRISPRs. This is very much a threshold that needs to be tested.
    next if($score<0.7);
    my $seqname=$$crispr[0][0];
    $score=sprintf("%.02f",$score);
    # CRISPR parent element
    print GFF join("\t",$seqname,"CG-Pipeline","ncRNA",$crisprStart,$crisprStop,$score,'+','.',"ID=$crisprId;Parent=$seqname;Name=$crisprId")."\n";

    for(my $i=0;$i<$numDrs;$i++){
      my($seqname,$start,$stop)=@{$$crispr[$i]};
      my($drStart,$drStop)=($start+1,$stop+1);
      my $sequence=substr($$contigs{$seqname},$start,$stop-$start+1);
      # Direct Repeat
      my $DRid=sprintf("DR%s",$drCounter++);
      print GFF join("\t",$seqname,"CG-Pipeline","direct_repeat",$drStart,$drStop,$score,'+','.',"ID=$DRid;Name=$DRid")."\n";
      # add a spacer between DRs
      if($i<$numDrs-1){
        my $nextDrStart=$$crispr[$i+1][1]+1;
        my($spacerStart,$spacerStop)=($drStop+1,$nextDrStart-1);
        my $spacerSeq=substr($$contigs{$seqname},$spacerStart-1,$spacerStop-$spacerStart+1);
        # Spacer
        my $spacerid=sprintf("spacer%s",$spacerCounter++);
        print GFF join("\t",$seqname,"CG-Pipeline","ncRNA",$spacerStart,$spacerStop,$score,'+','.',"ID=$spacerid;Name=$spacerid;Parent=$seqname")."\n";
      }
    }
    $crisprCounter++;
  }
  close GFF;
  logmsg "Done making GFF. Found ".($crisprCounter-1)." CRISPRs";
  return $gff;
}

# find the same small repeats throughout the genome to pick up more CIRSPRs
# TODO Do I need to find small repeats on reverse strand also?
sub findSequencePositionsInGenome{
  my($drSeqs,$contigs,$settings)=@_;
  #print Dumper $drSeqs;
  my %drPos;
  for my $contigid(keys(%$contigs)){
    my $contig=$$contigs{$contigid};
    for my $drSeq(@$drSeqs){
      my $rcDrSeq=$drSeq;
      $rcDrSeq=~tr/ATCGatcg/TAGCtagc/;
      $rcDrSeq=reverse($rcDrSeq);
      my @drPos;
      while($contig=~m/($drSeq)/g){
        my($start,$stop)=($-[1],$+[1]);
        push(@drPos,[$contigid,$start,$stop,1]);
      }
      #while($contig=~m/($rcDrSeq)/g){
      #  my($start,$stop)=($-[1],$+[1]);
      #  push(@drPos,[$contigid,$start,$stop,-1]);
      #}
      push(@{ $drPos{$drSeq} },@drPos);
    }
  }
  my $exhaustiveDrPos=pickUpMoreDrs(\%drPos,$contigs,$settings);
  return $exhaustiveDrPos;
  return \%drPos;
}

sub pickUpMoreDrs{
  my($drPos,$contigs,$settings)=@_;

  logmsg "Using blast to refine DR hits";

  # put the contigs back to a file and make a blast database
  my $db="$$settings{tempdir}/blastdb.fasta";
  AKUtils::printSeqsToFile($contigs,$db);
  command("legacy_blast.pl formatdb -i $db -p F 2>/dev/null",{stdout=>1,quiet=>1});

  # blast with each DR
  my $drSeqArr=filterRepeatSeqs([keys(%$drPos)],$settings);
  my %newDrPos;
  for my $drSeq(@$drSeqArr){
    my @newDrPos;
    my $blsResults=command("echo '$drSeq'|legacy_blast.pl blastall -p blastn -d $db -e 0.05 -W 7 -F F -m 8 -a $$settings{numcpus}",{stdout=>1,quiet=>1});
    my @blsResults=split(/\n/,$blsResults);
    logmsg "".scalar(@blsResults)."\t$drSeq";
    for(@blsResults){
      my($query,$hit,$identity,$length,$mismatch,$gap,$qStart,$qStop,$sStart,$sStop,$evalue,$score)=split(/\t/,$_);
      next if($length<$$settings{minCrisprSize} || $gap>0 || $mismatch>1 || $sStart>$sStop);
      push(@newDrPos,[$hit,$sStart-1,$sStop-1]);
    }
    $newDrPos{$drSeq}=\@newDrPos;
  }
  logmsg "Done";
  return \%newDrPos;
}

sub findSmallTandemRepeats{
  my($seqs,$settings)=@_;
  my $minCrisprSize=$$settings{minCrisprSize}||20;

  # TODO repeat over a few different window sizes
  my $Q=Thread::Queue->new;
  my $outQ=Thread::Queue->new;
  my @thr;
  $thr[$_]=threads->create(\&findRepeatsWorker,$Q,$outQ,$minCrisprSize,$settings) for(0..$$settings{numcpus}-1);
  
  my %drSeq; # direct repeat sequence possibilities
  while(my($seqid,$seq)=each(%$seqs)){
    my $contigLength=length($seq);
    my @window=();
    for(my $i=0;$i<$contigLength;$i+=$$settings{stepSize}){
      #next if($i<2900000 || $i>2940000);
      push(@window,[$i,$seqid,$seq]);
      while($Q->pending>($$settings{numcpus}*1000)){
        logmsg "Sleeping";
        sleep 1;
      }
      #logmsg "".$Q->pending." ".$outQ->pending if($i%100==0);
      if(@window>$$settings{numcpus}*10){ # enqueue if it is 10x CPUs
        $Q->enqueue(@window);
        @window=();
        #logmsg "Enqueued. ".$Q->pending." ".$outQ->pending;
      }
    }
    $Q->enqueue(@window);
    @window=();
  }
  # join and terminate
  $Q->enqueue(undef) for(@thr);
  $_->join for(@thr);

  my @tmp=($outQ->pending)?$outQ->extract(0,$outQ->pending):();
  $drSeq{$_}=1 for(@tmp);

  my $smallTandemRepeatSeqs=filterRepeatSeqs([keys(%drSeq)],$settings);
  return $smallTandemRepeatSeqs;
}

sub findRepeatsWorker{
  my($Q,$outQ,$minCrisprSize,$settings)=@_;
  my $TID="TID".threads->tid;
  my $tmpfifo="$$settings{tempdir}/fifo.$TID.tmp";
  mkfifo($tmpfifo,'0500') if(!-e $tmpfifo);
  my %seen;
  while(defined(my $tmp=$Q->dequeue)){
    my($i,$seqid,$seq)=@$tmp;
    my $start=$i;
    my $stop=$start+$$settings{windowSize}-1;
    logmsg "Looking for repeats at bp $i in $seqid" if($i%(1000*$$settings{stepSize})==0);

    my $subAsm=substr($seq,$start,$stop-$start+1);
    my $tmpRepeats=findRepeats($subAsm,$minCrisprSize,$settings);
    next if(!keys(%$tmpRepeats));
    for(keys(%$tmpRepeats)){
      next if($seen{$_}++);
      $outQ->enqueue($_);
    }
  }
  return 1;
}

sub findRepeats{
  my($subAsm,$minCrisprSize,$settings)=@_;
  my $TID="TID".threads->tid;
  my $tmpfifo="$$settings{tempdir}/fifo.$TID.tmp";
  my $tmpout=command("echo -e \">seq\n$subAsm\" > $tmpfifo & $$settings{repeatMatch} -E -t -n $minCrisprSize $tmpfifo 2>/dev/null",{quiet=>1,stdout=>1});
  
  # read the repeats file
  my %repeat;
  my $i=0;
  while($tmpout=~/\s*(.+)\s*\n/g){ # go through lines and trim
    next if($i++<2); # skip headers
    my($start1,$start2,$length)=split /\s+/, $1;
    next if($length>$$settings{maxCrisprSize});
    if($start2=~/(\d+)r/){ # revcom
      $start2=$1-1; # base 1 to base 0
      push(@{ $repeat{$start1} },[$start2,$length*-1]);
    } else { # not revcom
      push(@{ $repeat{$start1} },[$start2,$length]);
    }
  }

  # convert coordinates to sequences
  my %drSeq;
  while(my($start1,$repeatList)=each(%repeat)){
    for(@$repeatList){
      my($start,$length)=@$_;
      my $subseq;
      if($length<0){
        my $tmp=$start;
        $start=$start-abs($length)+1;
        $subseq=substr($subAsm,$start,abs($length));
        $subseq=~tr/ATCGatcg/TAGCtagc/;
        $subseq=reverse($subseq);
      } else {
        $subseq=substr($subAsm,$start,abs($length));
      }
      $drSeq{$subseq}++;
    }
  }
  #die "\n".Dumper keys(%drSeq) if(keys(%repeat)>2);

  return \%drSeq;
}

sub drAlignmentScore{
  my($crispr,$contigs,$settings)=@_;
  
  my $TID="TID".threads->tid;
  my $seqIn ="$$settings{tempdir}/in.$TID.fna";
  my $seqOut="$$settings{tempdir}/out.$TID.fna";

  # make input for alignment
  my $i=1;
  open(FNA,">",$seqIn) or die "Could not open $seqIn:$!";
  for(@$crispr){
    # get all start/stop from @$crispr
    my($seqname,$start,$stop)=@$_;
    # get substr/subseq
    my $sequence=substr($$contigs{$seqname},$start,$stop-$start+1);
    print FNA ">seq$i\n$sequence\n";
    $i++;
  }
  close FNA;

  # make multiple sequence alignment
  my $muscle=AKUtils::fullPathToExec("muscle");
  command("$muscle -in $seqIn -out $seqOut -quiet",{quiet=>1});

  # create an alignment score
  my $aln=Bio::AlignIO->new(-file=>$seqOut)->next_aln;
  my $identity=$aln->percentage_identity/100;
  my $alnScore=sprintf("%0.2f",1-((1-$identity)*2));
  return $alnScore;
}

# run a command
# settings params: warn_on_error (boolean)
#                  stdout to return stdout instead of error code
#                  quiet to not print the command unless an error/warning
# returns error code
# TODO run each given command in a new thread
sub command{
  my ($command,$settings)=@_;
  my $return;
  my $TID="TID".threads->tid;
  my $refType=ref($command);
  if($refType eq "ARRAY"){
    my $err=0;
    for my $c(@$command){
      $err+=command($c,$settings);
    }
    return $err;
  }
  logmsg "$TID RUNNING COMMAND\n  $command" if(!$$settings{quiet});
  if($$settings{stdout}){
    $return=`$command`;
  } else {
    system($command);
    $return=$?;
  }
  if($?){
    my $msg="ERROR running command $command\n  With error code $?. Reason: $!\n  in subroutine ".(caller(1))[3];
    logmsg $msg if($$settings{warn_on_error});
    die $msg if(!$$settings{warn_on_error});
  }
  return $return;
}


sub usage{
  my($settings)=@_;
  local $0=fileparse $0;
  "Finds CRISPRs in a fasta file
  Usage: $0 assembly.fasta -o outfile.sql
    -w window size in bp. Default: $$settings{windowSize}
    -s step size in bp. Default: $$settings{stepSize}
  "
}
