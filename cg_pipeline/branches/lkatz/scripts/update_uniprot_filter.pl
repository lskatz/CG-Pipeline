#!/usr/bin/env perl

#TODO pay attention to why things are not getting filtered
# The mechanism is there but the ID doesn't match up with locus name...?

# Filter out genes from a fasta file that belong to Eukaryotes
# Or, retain genes that belong to prokies
# author: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use LWP::Simple;
use FindBin;
use lib "$FindBin::RealBin/../lib";
$ENV{PATH} = "$FindBin::RealBin:".$ENV{PATH};
use AKUtils qw(logmsg);
use Bio::Perl;

$0=fileparse($0);
exit(main());

sub main{
  my $settings={
    url=>"ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs",
    filterOut=>[qw(yeast human humchr rice arath calbican cdlist celegans dicty fly hoxlist mgdtosp pkinfam pombe scorpktx)],
  };
  GetOptions($settings,qw(url=s infile=s@ outfile=s help tempdir=s));
  $$settings{tempdir}||=AKUtils::mktempdir();

  my $infile=$$settings{infile} or die "Error: no infile(s) given\n".usage($settings);
  my $outfile=$$settings{outfile} or die "Error: no outfile given\n".usage($settings);
  my $url=$$settings{url} or die "Internal error with URL parameter";
  my @filterOut=@{ $$settings{filterOut} };
  my %filterOutIndex;
    $filterOutIndex{$_}=1 for(@filterOut);

  # read which proteins to filter out
  my @proteinList=getBlacklist($url,\@filterOut,$settings);

  # post-process
  for (@proteinList){
    s/^\s+|\s+$//g;
    $_="" if(/\s+/);
    next if(!/\(|;/);
    # ($locus)=>$locus
    s/^\(|\)$//g;
    # $locus; => $locus
    s/;$//;
  }
  
  my %proteinList;
  $proteinList{$_}=1 for (@proteinList);
  logmsg "Done with filtering and indexing proteins to filter out";

  # read the infile(s) and write the filtered entries to the outfile
  logmsg "Writing the output file";
  my $out=Bio::SeqIO->new(-file=>">$outfile");
  for my $file(@$infile){
    my $in=Bio::SeqIO->new(-file=>$file);
    while(my $seq=$in->next_seq){
      next if($proteinList{$seq->id});
      $out->write_seq($seq);
    }
  }

  logmsg "Outfile is in $outfile";
  return 0;
}

sub getBlacklist{
  my($url,$filterOut,$settings)=@_;
  my $urlContent=get($url);
  my (@blacklist);
  for my $fileInfo(split(/\n/,$urlContent)){
    my $file=(split(/\s+/,$fileInfo))[-1];
    my($search,$ext)=fileparse($file,qw(.txt));
    $search=~s/\d+$//; # remove trailing numbers from the base search
    my $is_blacklisted=scalar(grep(/$search/i,@$filterOut));
    next if(!$is_blacklisted);

    #next if($search ne 'yeast'); # DEBUG

    my @blacklistGenes=genesListedInUniptrotTxt($file,$url,$settings);
    die "Did not blacklist any genes from $file" if (!@blacklistGenes);
    logmsg "$file yielded ".scalar(@blacklistGenes)." genes";
    push(@blacklist,@blacklistGenes);
  }
  return @blacklist;
}

# parse a text file on uniprot to find which genes are listed
sub genesListedInUniptrotTxt{
  my($file,$url,$settings)=@_;
  
  $file="$url/$file";
  my @locus;
  my $seen_header=0;
  for my $line(split(/\n/,get($file))){
    # test the first few fields
    for my $locus((split(/\s{2,}/,$line))[0..6]){
      if(is_locus($locus)){
        push(@locus,split(/\s*;\s*/,$locus));
      }
    }
  }
  return @locus;
}

# see if a locus matches some patterns
sub is_locus{
  my($locus)=@_;

  return 0 if(!defined($locus));
  return 0 if($locus=~/^\s*$/);
  # a locus does not start with a dash or underscore
  return 0 if($locus=~/^\s*[\-_]/);
  # a locus does not have lowercase, but if it does then it will have a digit
  # but limit the string length because there might be a huge chemical name
  #     e.g. 6-methopentyl.....5-N-kinase
  return 0 if($locus=~/[a-z]/ && $locus!~/.+\d/ && length($locus) > 8);

  return 1;
}

sub usage{
  my($settings)=@_;
  "Reads the uniprot/sprot/trembl databases and attempts to shrink
  them by including only prokaryotic genes
  Usage: $0 -i database.fasta [-i database2.fasta] -o outfile.fasta
  -i database.fasta
    The fasta file from uniprot, etc
  -o outfile
    The output fasta
  -u URL
    The basename URL for where genes are listed.
    Defaut: $$settings{url}
  -h
    This help menu
  ";
}

