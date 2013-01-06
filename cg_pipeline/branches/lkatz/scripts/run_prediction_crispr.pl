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

local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
exit(main());

sub main{
  local $0=fileparse $0;
  my $settings={
    appname=>"cgpipeline",
    minCrisprSize=>20,
    maxCrisprSize=>50,
  };
  GetOptions($settings,qw(outfile=s help windowSize=i stepSize=i tempdir=s minCrisprSize=i maxCrisprSize=i));
  die usage($settings) if(@ARGV<1 || $$settings{help});
  my $outfile=$$settings{outfile} || "$0.gff";
  my $infile=$ARGV[0];
  die "ERROR: Could not find infile $infile" if(!-e $infile);
  $$settings{tempdir}||=AKUtils::mktempdir();
  logmsg "Temporary directory is $$settings{tempdir}";
  logmsg "Detecting CRISPRs in $infile";

  my $contigs=AKUtils::readMfa($infile,{first_word_only=>1});
  my $smallTandemRepeatSeqs=findSmallTandemRepeats($contigs,$settings);
  my $positions=findSequencePositionsInGenome($smallTandemRepeatSeqs,$contigs,$settings);
  my $gff=positionsToGff($positions,$contigs,$settings);
  copy($gff,$outfile);
  logmsg "Finished. GFF is in $outfile";
  return 0;
}

sub positionsToGff{
  my($positions,$contigs,$settings)=@_;
  my $gff="$$settings{tempdir}/crispr.gff";
  open(GFF,">",$gff) or die "Could not open temporary GFF:$!";
  print GFF "##gff-version   3\n";
  # TODO print source line, describing the whole contig(s)
  my $crisprCounter=1;
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
      $score=$score*.9 if($spacerLength>2.5*$drLength || $spacerLength<0.6*$drLength);
    }

    # don't accept low-quality CRISPRs. This is very much a threshold that needs to be tested.
    next if($score<0.75);
    my $seqname=$$crispr[0][0];
    print GFF join("\t",$seqname,"CG-Pipeline","CRISPR",$crisprStart,$crisprStop,$score,'+','.',"ID=$crisprId")."\n";

    for(my $i=0;$i<$numDrs;$i++){
      my($seqname,$start,$stop)=@{$$crispr[$i]};
      my($drStart,$drStop)=($start+1,$stop+1);
      my $sequence=substr($$contigs{$seqname},$start,$stop-$start+1);
      print GFF join("\t",$seqname,"CG-Pipeline","DirectRepeat",$drStart,$drStop,$score,'+','.',"Parent=$crisprId;DR=$sequence")."\n";
      # add a spacer between DRs
      if($i<$numDrs-1){
        my $nextDrStart=$$crispr[$i+1][1]+1;
        my($spacerStart,$spacerStop)=($drStop+1,$nextDrStart-1);
        my $spacerSeq=substr($$contigs{$seqname},$spacerStart-1,$spacerStop-$spacerStart+1);
        print GFF join("\t",$seqname,"CG-Pipeline","Spacer",$spacerStart,$spacerStop,$score,'+','.',"Parent=$crisprId;Spacer=$spacerSeq")."\n";
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
  return \%drPos;
}

sub findSmallTandemRepeats{
  my($seqs,$settings)=@_;
  $$settings{windowSize}||=$$settings{maxCrisprSize}*5;
  $$settings{stepSize}||=$$settings{windowSize}-2*$$settings{maxCrisprSize};
  my $minCrisprSize=$$settings{minCrisprSize}||20;
  
  my %drSeq; # direct repeat sequence possibilities
  while(my($seqid,$seq)=each(%$seqs)){
    my $contigLength=length($seq);
    for(my $i=0;$i<$contigLength;$i+=$$settings{stepSize}){
    #for(my $i=2920000;$i<$contigLength;$i+=$$settings{stepSize}){
      my $start=$i;
      my $stop=$start+$$settings{windowSize}-1;
      logmsg "Finished with $i bp" if($i%100000==0);
      #if($i>2980000){
      #  logmsg "DEBUG";
      #  last;
      #}

      my $tmpIn=createSubassembly($seq,$start,$stop,$settings);
      my $tmp=findRepeats($tmpIn,$minCrisprSize,$settings);
      %drSeq=(%drSeq,%$tmp);
    }
  }

  return [keys(%drSeq)];
}

sub findRepeats{
  my($file,$minCrisprSize,$settings)=@_;
  my $tmpout="$$settings{tempdir}/repeat-match.out";
  my $repeatMatch=AKUtils::fullPathToExec("repeat-match");
  command("$repeatMatch -E -t -n $minCrisprSize '$file' 2>/dev/null | tail -n +3 | sort -k1n -k2n > $tmpout",{quiet=>1});
  #if(`wc -l $tmpout | awk '{print \$1}'`>5){system("cat $tmpout");die;}
  return {} if(`wc -l $tmpout | awk '{print \$1}'`<1);
  my $tmp=AKUtils::readMfa($file);
  my(undef,$seq)=each(%$tmp);
  
  # read the repeats file
  my %repeat;
  open(REPEAT,"<",$tmpout) or die "Internal error: $!";
  while(<REPEAT>){
    s/^\s+|\s+$//g; # trim
    my($start1,$start2,$length)=split /\s+/;
    next if($length>$$settings{maxCrisprSize});
    if($start2=~/(\d+)r/){ # revcom
      $start2=$1-1; # base 1 to base 0
      push(@{ $repeat{$start1} },[$start2,$length*-1]);
    } else { # not revcom
      push(@{ $repeat{$start1} },[$start2,$length]);
    }
  }

  # add the start1 site as a repeat too
  while(my($start1,$repeatList)=each(%repeat)){
    my $length=abs($$repeatList[0][1]);
    #push(@{ $repeat{$start1} },[$start1,$length]);
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
        $subseq=substr($seq,$start,abs($length));
        $subseq=~tr/ATCGatcg/TAGCtagc/;
        $subseq=reverse($subseq);
      } else {
        $subseq=substr($seq,$start,abs($length));
      }
      $drSeq{$subseq}++;
    }
  }
  #die "\n".Dumper keys(%drSeq) if(keys(%repeat)>2);

  return \%drSeq;
}

sub drAlignmentScore{
  my($crispr,$contigs,$settings)=@_;
  
  my $seqIn ="$$settings{tempdir}/in.fna";
  my $seqOut="$$settings{tempdir}/out.fna";

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

sub createSubassembly{
  my($seq,$start,$stop,$settings)=@_;
  my $subasm="$$settings{tempdir}/asm.$start.$stop.fasta";
  open(ASM,">",$subasm) or die "Could not open sub assembly for writing: $!";
  print ASM join("_",">subasm",$start,$stop)."\n".substr($seq,$start,($stop-$start+1))."\n";
  close ASM;
  return $subasm;
}

# run a command
# settings params: warn_on_error (boolean)
# returns error code
# TODO run each given command in a new thread
sub command{
  my ($command,$settings)=@_;
  my $refType=ref($command);
  if($refType eq "ARRAY"){
    my $err=0;
    for my $c(@$command){
      $err+=command($c,$settings);
    }
    return $err;
  }
  logmsg "RUNNING COMMAND\n  $command" if(!$$settings{quiet});
  system($command);
  if($?){
    my $msg="ERROR running command $command\n  With error code $?. Reason: $!\n  in subroutine ".(caller(1))[3];
    logmsg $msg if($$settings{warn_on_error});
    die $msg if(!$$settings{warn_on_error});
  }
  return $?;
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
