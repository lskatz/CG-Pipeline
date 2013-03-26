#!/usr/bin/env perl

# Predict CRISPRs in a fasta file.
# Author: Lee Katz <lkatz@cdc.gov>

# Note: if any other CRISPR detection programs are added to CG-Pipeline,
# they should be added here.

use strict;
use warnings;
use FindBin;
use lib "$FindBin::RealBin/../lib";
$ENV{PATH} = "$FindBin::RealBin:".$ENV{PATH};
use AKUtils;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use File::Copy qw/copy move/;
use Bio::AlignIO;
use List::Util qw/sum min max/;

my $directRepeatId=0; # this is a global var to keep track of unique IDs

$0=fileparse $0;
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
exit(main());

sub logmsg{print STDERR "$0: ".(caller(1))[3].": @_\n";}
sub main{
  my $settings={
    appname=>"cgpipeline",
  };
  $settings=AKUtils::loadConfig($settings);
  GetOptions($settings,qw(help tempdir=s crtJar=s));
  $$settings{tempdir}||=AKUtils::mktempdir($settings);
  $$settings{crtJar}||=$$settings{prediction_crt_jar};

  #$$settings{pilerCr}=AKUtils::fullPathToExec("pilercr");
  my $outfile=$$settings{outfile} || "$0.gff";
  die usage() if(@ARGV<1);
  my @assembly=@ARGV;
  
  my @gff;
  for my $assembly(@assembly){
    logmsg "$assembly";
    my $crtGff=runCrt($assembly,$settings);
    #my $pilerCrOut=runPilerCr($assembly,$settings);
    #my $gff=pilerCrToGff($pilerCrOut,$settings);

    my $gff=$crtGff; # combine the gffs somehow
    push(@gff,$gff);
  }

  # write the composite GFF, sorting by contig and then start/stop
  print "##gff-version  3\n";
  open(GFF,"sort -k1,1 -k4,4n -k5,5nr ".join(" ",@gff)." |") or die "Could not open GFFs:$!";
  while(<GFF>){
    print;
  }
  close GFF;
  
  return 0;
}
sub runCrt{
  my($assembly,$settings)=@_;
  my $crtJar=$$settings{crtJar}||"";
  die "ERROR Cannot find CRT jar file at $crtJar. If not set, then add a value for prediction_crt_jar to your config file" if(!-f $crtJar);
  my $basename=fileparse $assembly;
  my $crtOut="$$settings{tempdir}/$basename.crt.out";
  my $java=AKUtils::fullPathToExec("java");
  my $seqs=AKUtils::readMfa($assembly,{first_word_only=>0});
  # CRT has a problem with multifasta, so run it once per contig
  mkdir "$$settings{tempdir}/gff";
  while(my($id,$sequence)=each(%$seqs)){
    my $contig="$$settings{tempdir}/tmpcontig.fasta";
    AKUtils::printSeqsToFile({$id=>$sequence},$contig);
    my $command="$java -cp $crtJar crt -minNR 2 -maxSL 100 '$contig' '$crtOut' 1>/dev/null 2>&1";
    system($command);
    die if $?;
    my $crtGff=crtToGff($crtOut,$settings);
    move("$crtGff","$$settings{tempdir}/gff/$id.gff"); die $! if $?;
  }
  system("cat $$settings{tempdir}/gff/*.gff > $$settings{tempdir}/$basename.gff"); die if $?;
  return "$$settings{tempdir}/$basename.gff";
}

sub crtToGff{
  my($crtOut,$settings)=@_;
  my $gff="$crtOut.gff";
  open(CRT,"<",$crtOut) or die "Could not read $crtOut:$!";
  open(GFF,">",$gff) or die "Could not write to $gff:$!";
  my $drId=$directRepeatId; # value is from a global scope
  my $seqname="UNDEFINED";
  while(<CRT>){
    if(/ORGANISM:\s+(.+)\s*$/){
      $seqname=$1;
    }
    if(/CRISPR (\d+)/){
      my $crisprId=$1;
      my $crisprStart=0;
      my $crisprStop=0;
      while(<CRT>){
        next if(/^POSITION/);
        next if(/^\-\-/);
        last if(/^Repeats/);
        chomp;
        my($drStart,$DR,$spacer,$lengths)=split /\t+/;
        my $drLength=length($DR);
        my $drStop=$drStart+$drLength-1;
        $crisprStart||=$drStart;
        if($spacer){
          my $spacerLength=length($spacer);
          my $spacerStart=$drStop+1;
          my $spacerStop=$spacerStart+$spacerLength-1;
        } else {
          $crisprStop=$drStop;
        }
        $drId++;
        print GFF join("\t",$seqname,"CrisprRecognitionTool","direct_repeat",$drStart,$drStop,'.','.','.',"Parent=CRISPR$crisprId;Name=DR$drId;ID=DR$drId")."\n";
      }
      next if(!$crisprStart);
      print GFF join("\t",$seqname,"CrisprRecognitionTool","repeat_region",$crisprStart,$crisprStop,'.','.','.',"Name=CRISPR$crisprId;Id=CRISPR$crisprId")."\n";
    }
  }
  $directRepeatId=$drId; # return the value back to global scope
  close CRT;
  close GFF;
  return $gff;
}

sub runPilerCr{
  my($assembly,$settings)=@_;

  my $basename=fileparse $assembly;
  my $pilerCrOut="$$settings{tempdir}/$basename.pilerCr.out";
  my $pilerCrDr="$$settings{tempdir}/$basename.pilerCr.dr.fasta";
  # pilercr -in P_aeruginoasP14.fasta -out /dev/stdout -seq pilecr.seq -noinfo
  system("$$settings{pilerCr} -in '$assembly' -out '$pilerCrOut' -seq '$pilerCrDr' -noinfo -quiet");
  die if $?;

  return ($pilerCrOut,$pilerCrDr) if wantarray;
  return $pilerCrOut;
}

sub pilerCrToGff{
  my($pilerCrOut,$settings)=@_;
  my $gff="$pilerCrOut.gff";
  open(PILERCR,"<",$pilerCrOut) or die "Could not open piler-cr file $pilerCrOut:$!";
  open(GFF,">",$gff) or die "Could not open $gff for writing:$!";
  my $stage="";
  while(my $line=<PILERCR>){
    if($line=~/DETAIL REPORT/){
      $stage="detail";
    }
    elsif($line=~/SUMMARY BY SIMILARITY/){
      $stage="similarity";
    }
    elsif($line=~/SUMMARY BY POSITION/){
       $stage="position"; 
    }
    # detailed DRs
    elsif($stage eq "detail" && $line=~/Array\s+(\d+)/){
      my $crisprId=$1;
      my $seqname=<PILERCR>; chomp($seqname);
      $seqname=~s/\s+.*//; # remove anything after whitespace
      $seqname=~s/^>//; # remove > in defline
      my $undef=<PILERCR> for(1..3); # discard: blank, headers, ===
      my($start,$stop,@score);
      while(<PILERCR>){
        last if /===/;
        my @F=split /\s+/; @F=grep(!/^\s*$/,@F);
        my $drStop=($F[0]+$F[1]-1);
        print GFF join("\t",$seqname,"Piler-CR","direct_repeat",$F[0],$drStop,$F[2],'.','.',"Parent=CRISPR${crisprId}")."\n";

        # CRISPR region data
        push(@score,$F[2]);
        $start=$F[0] if(!defined($start) || $start>$F[0]);
        $stop=$drStop if(!defined($stop) || $stop<$drStop);
      }
      my $score=sum(@score)/@score;
      print GFF join("\t",$seqname,"Piler-CR","repeat_region",$start,$stop,$score,'.','.',"ID=CRISPR${crisprId};Name=CRISPR${crisprId}")."\n";
    }

  }
  close GFF;
  close PILERCR;
  return $gff;
}

sub usage{
  my($settings)=@_;
  local $0=fileparse $0;
  "Finds CRISPRs in a fasta file
  Usage: $0 assembly.fasta > outfile.gff
    -t tempdir [optional]
  "
}
