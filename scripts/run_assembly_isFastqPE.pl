#!/usr/bin/env perl
# Sees if a fastq file is Paired End or not

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/basename/;

$0=basename $0;
sub logmsg{print STDERR "@_\n";}
exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(checkFirst=i help wantformat)) or die "$!\n".usage();
  die usage() if($$settings{help} || !@ARGV);

  # if checkFirst is undef or 0, this will cause it to check at least the first 20 entries.
  # 200 reads is probably enough to make sure that it's shuffled (1/2^100 chance I'm wrong)
  $$settings{checkFirst}||=200;
  $$settings{checkFirst}=200 if($$settings{checkFirst}<2);

  my $file=shift @ARGV;

  my ($is_PE,$format)=is_fastqInterlevedPE($file,$settings);
  print "$is_PE\n";
  print "$format\n" if($$settings{wantformat});
  return 0;
}


# See whether a fastq file is paired end or not. It must be in a velvet-style shuffled format.
# In other words, the left and right sides of a pair follow each other in the file.
# params: fastq file and settings
# fastq file can be gzip'd
# settings:  checkFirst is an integer to check the first X deflines
sub is_fastqInterlevedPE($;$){
  my($fastq,$settings)=@_;

  # get the deflines
  my @defline;
  my $numEntries=0;
  my $i=0;
  my $fp;
  if($fastq=~/\.gz$/){
    open($fp,"gunzip -c '$fastq' |") or die "Could not open $fastq for reading: $!";
  }else{
    open($fp,"<",$fastq) or die "Could not open $fastq for reading: $!";
  }
  my $discard;
  while(my $defline=<$fp>){
    next if($i++ % 4 != 0);
    chomp($defline);
    $defline=~s/^@//;
    push(@defline,$defline);
    $numEntries++;
    last if($numEntries > $$settings{checkFirst});
  }
  close $fp;

  # it is paired end if it validates with any naming system
  my $is_pairedEnd=0;
  my $format="NotPE";
  if(_is_fastqPESra(\@defline,$settings)){
    $is_pairedEnd=1;
    $format="SRA";
  } elsif (_is_fastqPECasava18(\@defline,$settings)){
    $is_pairedEnd=1;
    $format="Casava1.8";
  } elsif (_is_fastqPECasava17(\@defline,$settings)){
    $is_pairedEnd=1;
    $format="Casava1.7";
  }
  return ($is_pairedEnd,$format) if wantarray;
  return $is_pairedEnd;
}

sub _is_fastqPESra{
  my($defline,$settings)=@_;
  my @defline=@$defline; # don't overwrite $defline by mistake

  for(my $i=0;$i<@defline-1;$i+=2){
    my($genome,$info1,$info2)=split(/\s+/,$defline[$i]);
    if(!$info2){
      return 0;
    }
    my($instrument,$flowcellid,$lane,$x,$y,$X,$Y)=split(/:/,$info1);
    my($genome2,$info3,$info4)=split(/\s+/,$defline[$i+1]);
    my($instrument2,$flowcellid2,$lane2,$x2,$y2,$X2,$Y2)=split(/:/,$info3);
    $_||="" for($X,$Y,$X2,$Y2); # these variables might not be present
    if($instrument ne $instrument2 || $flowcellid ne $flowcellid2 || $lane ne $lane2 || $x ne $x2 || $y ne $y2 || $X ne $X2 || $Y ne $Y2){
      return 0;
    }
  }
  return 1;
}
sub _is_fastqPECasava18{
  my($defline,$settings)=@_;
  my @defline=@$defline;

  for(my $i=0;$i<@defline-1;$i+=2){
    my($instrument,$runid,$flowcellid,$lane,$tile,$x,$yandmember,$is_failedRead,$controlBits,$indexSequence)=split(/:/,$defline[$i]);
    $yandmember="" if(!defined($yandmember));
    my($y,$member)=split(/\s+/,$yandmember);

    my($inst2,$runid2,$fcid2,$lane2,$tile2,$x2,$yandmember2,$is_failedRead2,$controlBits2,$indexSequence2)=split(/:/,$defline[$i+1]);
    $yandmember2="" if(!defined($yandmember2));
    my($y2,$member2)=split(/\s+/,$yandmember2);

    # If there isn't even a member, then this is not the correct format.
    return 0 if(!$member && !$member2);

    # Instrument, etc must be the same.
    # The member should be different, usually "1" and "2"
    $_||="" for($instrument,$inst2,$runid,$runid2,$flowcellid,$fcid2,$tile,$tile2);
    $_||=-1 for($member,$member2);
    if($instrument ne $inst2 || $runid ne $runid2 || $flowcellid ne $fcid2 || $tile ne $tile2 || $member>=$member2){
      return 0;
    }
  }
  return 1;
}

# This format is basically whether the ends of the defline alternate 1 and 2.
sub _is_fastqPECasava17{
  my($defline,$settings)=@_;
  my @defline=@$defline;
  for(my $i=0;$i<@defline-1;$i+=2){
    # Get each member number but return false if it doesn't even exist.
    my ($member1,$member2);
    if($defline[$i] =~ m/(\d+)$/){
      $member1=$1;
    } else {
      return 0;
    }
    if($defline[$i+1] =~ /(\d+)$/){
      $member2=$1;
    } else {
      return 0;
    }

    # The test is whether member1 is less than member2.
    # They can't be equal either.
    if($member1 >= $member2){
      return 0;
    }

    # Currently I don't believe that the members should be double digits which I guess could change
    return 0 if($member2 > 9 || $member1 > 9);
  }

  return 1;
}

sub usage{
  "$0: checks to see if a reads file is paired-end or not. Outputs 1 for yes or 0 for no.
  Usage: $0 file.fastq[.gz]
  --checkFirst [200]  How many deflines to check to make sure it is PE
  --wantformat        Use this flag to also print which defline format the file is
  "
}

