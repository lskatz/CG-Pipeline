#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Perl;
use Data::Dumper;
use LWP::Simple;
use threads;
use Thread::Queue;
# TODO multithread SeqIO/in
# TODO buffer the output
# TODO download the correct genus files
sub logmsg{print STDERR "@_\n";}
exit main();
sub main{
  my $settings={};
  GetOptions($settings,qw(help tempdir=s)) or die;
  die usage() if($$settings{help});
  die "ERROR: need fasta files\n".usage() if(!@ARGV);
  my @fasta=@ARGV;
  my @genera=getGenera($settings);
  for my $fasta(@fasta){
    printNewFasta($fasta,\@genera,$settings);
  }
  return 0;
}
sub getGenera{
  my($settings)=@_;
  # download the taxonomy html files
  my $html;
  for my $url(qw(http://www.bacterio.net/-ac.html http://www.bacterio.net/-dl.html http://www.bacterio.net/-mr.html http://www.bacterio.net/-sz.html)){
    logmsg "Downloading $url";
    $html.=get($url);
  }
  logmsg "Finding genus names";
  my @genus;
  for my $line(split(/\n/,$html)){
    while($line=~/'genusspecies'>(\w+)</gi){
      push(@genus,$1);
    }
  }
  return @genus;
}
sub printNewFasta{
  my($fastain,$genus,$settings)=@_;
  logmsg "Reading $fastain and writing to stdout";
  my $out=Bio::SeqIO->new(-format=>"fasta");
  my $in=Bio::SeqIO->new(-file=>$fastain);
  my $genusRegex=join("|",@$genus);
  my $i=0;
  while(my $seq=$in->next_seq){
    logmsg $i." entries read." if(++$i % 10000 == 0);
    my $id=$seq->id." ".$seq->desc;
    if($id=~/$genusRegex/){
      my($seqid,$product,$gene)=("","","");
      if($id=~/^\s*(\S+)\s+(.*?)\s*\w+=/){
        ($seqid,$product)=($1,$2);
      }
      if($id=~/GN=(\S+)/){
        $gene=$1;
      }
      $seq->id($seqid);
      $seq->desc("~~~$gene~~~$product");
      $out->write_seq($seq);
    }
  }
}
sub usage{
  "Filters a fasta file with the genera in the bacterial kingdom
  Usage: $0 uniprot.fasta > filtered.fasta
  "
}

