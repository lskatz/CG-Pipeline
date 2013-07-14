#!/usr/bin/env perl
# Mass-trims the 5' and 3' ends of a read set based on average nt percentages
# (in the future) and maybe based on other criteria

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
use Data::Dumper;

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH=*STDERR; print $FH "$0: ".(caller(1))[3].": @_\n";}

exit(main());

sub main() {
  GetOptions($settings, qw(tempdir=s cutoffReadNum=i readlength=i fastqqc-verbose badSlope=s));
  die usage() if(@ARGV<1);
  $$settings{badSlope}||=0.7; # this is kind of empirical for a cutoff
  $$settings{cutoffReadNum}=0 if(!defined $$settings{cutoffReadNum});
  $$settings{readlength}||=999;
  logmsg "Read length was given as $$settings{readlength}bp. NOTE: Reads longer than $$settings{readlength} will be trimmed regardless of quality";
  my @file=@ARGV;

  for my $file (@file){
    my $trimmedFile=trimUpAndDown($file,$settings);
  }

  return 0;
}

# Trim upstream and downstream based on the stability of ATCG percentages.
# Uses fastqqc internally
sub trimUpAndDown{
  my($file,$settings)=@_;

  my $numSlopesToBeGood=2;
  my $badSlope=$$settings{badSlope}; 

  my $fastqqc=AKUtils::fullPathToExec("fastqqc");
  my $fastqqcXopts="";
    $fastqqcXopts.="-c $$settings{cutoffReadNum} " if($$settings{cutoffReadNum});
    $fastqqcXopts.="-r $$settings{readlength} " if($$settings{readlength});
  logmsg "Running $fastqqc on $file and looking for balanced base percentages along $numSlopesToBeGood bases of the reads";
  # 5' and 3' trimming step
  my(@a,@t,@c,@g,@q);
  if(is_fastqGz($file,$settings)){
    open(STREAM,"gunzip -c '$file' | $fastqqc $fastqqcXopts /dev/stdin |") or die "Could not run fastqqc: $!";
  } else {
    open(STREAM,"$fastqqc $fastqqcXopts '$file' |") or die "Could not run fastqqc: $!";
  }
  while(<STREAM>){
    if(/^pos/){
      while(<STREAM>){
        last if(/Ns per read/);
        #print if($$settings{"fastqqc-verbose"});
        chomp;
        my($pos,$a,$c,$g,$t,$n,$q)=split /\t/; 
        push(@a,$a); push(@c,$c); push(@g,$g); push(@t,$t);
        push(@q,$q);
      }
    }
  } 
  close STREAM;
  if($$settings{"fastqqc-verbose"}){
    fastqqcOutput(\@a,\@c,\@g,\@t,\@q,$settings);
    return 0;
  }

  # find where to make the cut
  my ($upstreamTrim,$downstreamTrim);
  # 5' upstream
  my $numGoodInARow=0;
  for(my $i=0;$i<@a-1;$i++){
    my $g=sprintf("%.02f",abs($g[$i+1]-$g[$i])); 
    my $t=sprintf("%.02f",abs($t[$i+1]-$t[$i])); 
    my $a=sprintf("%.02f",abs($a[$i+1]-$a[$i])); 
    my $c=sprintf("%.02f",abs($c[$i+1]-$c[$i])); 
    my $slope=sprintf("%.02f",($a+$c+$t+$g)/4); 
    #print join("\t",$i,"$a[$i] $a","$c[$i] $c","$g[$i] $g","$t[$i] $t",$q[$i],$slope,"<=")."\n";
    # If the difference b/n this % and the next bp's % is large or if this bp has low qual, 
    # then it is not good.
    if($slope>$badSlope || $q[$i]==-64){ 
      $numGoodInARow=0;
    }else{
      $numGoodInARow++;
    } 
    if($numGoodInARow>=$numSlopesToBeGood){
      $upstreamTrim=$i-$numSlopesToBeGood;
      last;
    }
  }

  # 3' downstream
  $numGoodInARow=0;
  @a=reverse(@a);@c=reverse(@c);@t=reverse(@t);@g=reverse(@g);@q=reverse(@q);
  for(my $i=0;$i<@a-1;$i++){
    my $g=sprintf("%.02f",abs($g[$i+1]-$g[$i])); 
    my $t=sprintf("%.02f",abs($t[$i+1]-$t[$i])); 
    my $a=sprintf("%.02f",abs($a[$i+1]-$a[$i])); 
    my $c=sprintf("%.02f",abs($c[$i+1]-$c[$i])); 
    my $slope=sprintf("%.02f",($a+$c+$t+$g)/4); 

    if($slope>$badSlope || $q[$i]==-64){ # or a zero quality
      $numGoodInARow=0;
    }else{
      $numGoodInARow++;
    } 
    if($numGoodInARow>=$numSlopesToBeGood){
      $downstreamTrim=$i-$numSlopesToBeGood;
      last;
    }
  }
  close STREAM;

  my $start=$upstreamTrim; $start=0 if($start<0);
  # length to keep is the orig length minus the trimming lengths
  my $length=scalar(@a)-($downstreamTrim+$start);
  logmsg "Slopes were good starting at $upstreamTrim from 5' and $downstreamTrim from 3'. Therefore I will trim at $start and retain $length length reads.";
  if(is_fastqGz($file,$settings)){
    open(IN,"gunzip -c '$file' |") or die "Could not open $file for gunzipping/reading:$!";
  } else {
    open(IN,"<",$file) or die "Could not open $file for reading:$!";
  }
  my $i=0;
  while(<IN>){
    my $mod=$i%4;
    if($mod==1 || $mod==3){
      chomp;
      print substr($_,$start,$length)."\n";
    } else {
      print $_;
    }
    $i++;
  }
  close IN;
  return $i;
}

sub fastqqcOutput{
  my($a,$c,$g,$t,$q,$settings)=@_;
  my $numBases=@$a;
  print join("\t",qw(nt a c g t q slope))."\n";
  for(my $i=0;$i<$numBases;$i++){
    my $slope="-";
    if($i>0){
      $slope=(abs($$a[$i+1]-$$a[$i])+abs($$c[$i+1]-$$c[$i])+abs($$g[$i+1]-$$g[$i])+abs($$t[$i+1]-$$t[$i]))/4;
    }
    print join("\t",$i+1,$$a[$i],$$c[$i],$$g[$i],$$t[$i],$$q[$i],$slope)."\n";
  }
  return $numBases;
}

sub is_fastqGz($;$){
  my($file,$settings)=@_;
  return 1 if(ext($file)=~/^(fastq\.gz)$/i);
  return 0;
}
# determine file types
sub ext($){ # return the extension of a filename
  my($file)=@_;
  my($name,$path,$ext)=fileparse($file,qw(.fastq.gz .sff .fasta .fna .fa .fas .fastq .ffn .mfa .fq));
  if(!$ext){
    $ext=$file;
    $ext=~s/.+(\..*$)/\1/;
  }
  $ext=~s/^.//; # remove that dot
  return lc($ext);
}

sub usage{
  "Usage: $0 reads.fastq > reads.trimmed.fastq
    -c 9999 number of reads to sample before deciding on an average quality.
      0 for all (not recommended because it is very slow and the quality will not be too different)
    -r 999 the length of all reads. There is no effect of going over the read length, but it is bad if you go under
    --fastqqc-verbose print the fastqqc report and exit
    -b 0.7 the percentage change in average ATCG content that is considered to be a bad skew (default: 0.7%)
  "
}
