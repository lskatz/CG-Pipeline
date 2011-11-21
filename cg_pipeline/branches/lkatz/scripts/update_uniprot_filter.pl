#!/usr/bin/env perl

#TODO pay attention to why things are not getting filtered
# The mechanism is there but the ID doesn't match up with locus name...?

# Filter out genes from a fasta file that belong to Eukaryotes
# author: Lee Katz <lkatz@cdc.gov>

package PipelineRunner;
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
    appname=>"cgpipeline",
    url=>"ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs",
    filterOut=>[qw(yeast human humchr rice arath calbican cdlist celegans dicty fly hoxlist mgdtosp pkinfam pombe scorpktx)],
  };
  GetOptions($settings,qw(url=s infile=s@ outfile=s help tempdir=s));
  $$settings{tempdir}||=AKUtils::mktempdir();

  my $infile=$$settings{infile} or die "Error: no infile(s) given\n".usage($settings);
  my $outfile=$$settings{outfile} or die "Error: no outfile given\n".usage($settings);
  my $url=$$settings{url} or die "Internal error with URL parameter";
  my @filterOut=@{ $$settings{filterOut} };

  # read which proteins to filter out
  my @proteinList=getBlacklist($url,\@filterOut,$settings);

  my %proteinList;
  $proteinList{$_}=1 for (@proteinList);
  logmsg "Done with filtering and indexing proteins to filter out";

  # read the infile(s) and write the filtered entries to the outfile
  logmsg "Writing the output file";
  my $seqCount=0;
  my $out=Bio::SeqIO->new(-file=>">$outfile");
  for my $file(@$infile){
    my $in=Bio::SeqIO->new(-file=>$file);
    while(my $seq=$in->next_seq){
      next if($proteinList{$seq->id});
      $out->write_seq($seq);
      logmsg "Finished with $seqCount sequences" if(++$seqCount%1000000==0);
    }
  }

  logmsg "Outfile is in $outfile";
  return 0;
}

sub getBlacklist{
  my($url,$filterOut,$settings)=@_;
  my (@blacklist);
  my $tmpfile="$$settings{tempdir}/blacklist.txt";
  
  # if the blacklist already exists, just read the file to save time
  if(-e $tmpfile){
    logmsg "Blacklist already exists at $tmpfile; reading instead of recreating it";
    open(FILE,$tmpfile) or die "Could not open $tmpfile because $!";
    while(<FILE>){
      chomp;
      push(@blacklist,$_);
    }
    close FILE;
    return @blacklist;
  }

  my $urlContent=get($url);
  for my $fileInfo(split(/\n/,$urlContent)){
    my $file=(split(/\s+/,$fileInfo))[-1];
    my($search,$ext)=fileparse($file,qw(.txt));
    $search=~s/\d+$//; # remove trailing numbers from the base search
    my $is_blacklisted=scalar(grep(/$search/i,@$filterOut));
    next if(!$is_blacklisted);

    my @blacklistGenes=genesListedInUniptrotTxt($file,$url,$settings);
    warn "Did not blacklist any genes from $file" if (!@blacklistGenes);
    logmsg "$file yielded ".scalar(@blacklistGenes)." genes";
    push(@blacklist,@blacklistGenes);
  }
  # post-process
  for (@blacklist){
    s/^\s+|\s+$//g;
    $_="" if(/\s+/);
    next if(!/\(|;/);
    # ($locus)=>$locus
    s/^\(|\)$//g;
    # $locus; => $locus
    s/;$//;
  }
  @blacklist=grep(!/^\s*$/,@blacklist);

  # write the blacklist to file to cache it
  # /(sp|tr)\|[\dA-Z]{6}/
  open(FILE,">",$tmpfile) or die "Could not open $tmpfile for writing because $!";
  print FILE "$_\nsp|$_\ntr|$_\n" for(@blacklist);
  close FILE;
  logmsg "Cached the blacklist at $tmpfile";
  
  return @blacklist;
}

# download and parse a text file on uniprot to find which genes are listed
sub genesListedInUniptrotTxt{
  my($file,$url,$settings)=@_;
  
  $file="$url/$file";
  my @locus;
  for my $line(split(/\n/,get($file))){
    # test the first few fields
    for my $locus(split(/\s{2,}/,$line)){
      if(is_locus($locus)){
        my @theseLoci=split(/\s*;\s*/,$locus);
        for my $l(@theseLoci){
          push(@locus,$l) if(is_locus($l));
        }
      }
    }
  }
  return @locus;
}

# see if a locus matches some patterns
# Every locus is a 6-digit thing
#   grep ">" uniprot_sprot_trembl.fasta|perl -lane 'die "$_" if(!/(sp|tr)\|[\dA-Z]{6}/)'
sub is_locus{
  my($locus)=@_;
  return 0 if(!defined($locus));
  return 0 if($locus!~/[0-9A-Z]{6}/);
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
  -t directory
    a temporary directory. Default: a newly created dir under /tmp/
  -h
    This help menu
  ";
}

