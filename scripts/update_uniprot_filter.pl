#!/usr/bin/env perl

# Filter out genes from a fasta file that don't belong to bacteria
# author: Lee Katz <lkatz@cdc.gov>

package PipelineRunner;
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Fetch;
use Data::Dumper;
use LWP::Simple;
use FindBin;
use lib "$FindBin::RealBin/../lib";
$ENV{PATH} = "$FindBin::RealBin:".$ENV{PATH};
use AKUtils qw(logmsg);
use Bio::Perl;
use Data::Dumper;

$0=fileparse($0);
exit(main());

sub main{
  my $settings={
    appname=>"cgpipeline",
  };
  GetOptions($settings,qw(url=s@ infile=s@ outfile=s help tempdir=s));
  $$settings{tempdir}||=AKUtils::mktempdir();
  $$settings{url}=[qw(ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_bacteria.dat.gz ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_bacteria.dat.gz)] if(!$$settings{url});

  my $infile=$$settings{infile} or die "Error: no infile(s) given\n".usage($settings);
  my $outfile=$$settings{outfile} or die "Error: no outfile given\n".usage($settings);
  my $url=$$settings{url};

  logmsg "Downloading and uncompressing dat files from ".join(", ".@{$$settings{url}});
  my $file=downloadDats($url,$settings);
  logmsg "Parsing for the correct IDs to keep";
  my $bacterialIds=findBacterialIds($file,$settings);
  logmsg "Printing to file";
  my $numSeqs=printFilteredFile($bacterialIds,$infile,$outfile,$settings);

  logmsg "Done! $numSeqs sequences printed to $outfile";

  return 0;
}

sub downloadDats{
  my($url,$settings)=@_;
  my @file=();
  my $dir=".";
  for (@$url){
    my $path="$dir/".fileparse($_);
    my $datpath="$dir/".fileparse($_,'.gz');
    # download the dat file if it doesn't exist in gz or uncompressed form
    if(! -f $path && ! -f $datpath){
      logmsg "Downloading $_ => $path";
      system("wget '$_' -O '$path'"); die if $?;
    }
    if(! -f $datpath){
      logmsg "Uncompressing $path";
      system("gunzip '$path'"); die if $?;
    }
    push(@file,$datpath);
  }
  return \@file;
}

sub findBacterialIds{
  my($dat,$settings)=@_;
  my %id;
  my $i=0;
  for my $datfile (@$dat){
    # using grep because it is probably faster than going line by line and checking with perl regex
    open(DAT,"grep '^ID' '$datfile' |") or die "ERROR: could not open $datfile: $!";
    while(<DAT>){
      my $id=(split/\s+/,$_)[1];
      $id{$id}=1;
      $i++;
      if($i % 1000000 == 0){
        $|++;
        logmsg "Finished reading $i IDs from $datfile";
        $|--;
      }
    }
    close DAT;
  }
  return \%id;
}

sub printFilteredFile{
  my ($bacterialIds,$infile,$outfile,$settings)=@_;
  my $seqout=Bio::SeqIO->new(-file=>">$outfile");
  my $seqCount=0;
  my $i=0;
  for my $in(@$infile){
    my $seqin=Bio::SeqIO->new(-file=>$in);
    while(my $seq=$seqin->next_seq){
      my $id=(split(/\|/,$seq->id))[2];
      next if(!$$bacterialIds{$id});
      $seqout->write_seq($seq);
      $seqCount++;
      $i++;
      
      if($i % 1000000 == 0){
        logmsg "Searched through $i sequence records";
      }
    }
  }
  my $numFiltered=$i-$seqCount;
  logmsg "Done. Kept $seqCount sequences out of $i total; filtered $numFiltered sequences";
  return $seqCount;
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
    Defaut: ".join("\t",@{ $$settings{url} })."
  -h
    This help menu
  ";
}

