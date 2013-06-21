#!/usr/bin/env perl

# author: Jay Humphrey

use warnings;
use strict;
use File::Basename;
use FindBin;
use lib "$FindBin::RealBin/../lib";
use AKUtils;
use Bio::Seq;
use Bio::Seq::RichSeq;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Data::Dumper;
use Getopt::Long;

exit(&main());

sub main(){
  die usage() if(@ARGV<1);
  my $settings;
  $settings = AKUtils::loadConfig($settings);
  GetOptions($settings,('gb=s','fasta=s','faa=s')) or die;
  if(!defined($$settings{'gb'}) or !defined($$settings{'fasta'})){die usage();}
  printf STDERR "reading %s... ",$$settings{'gb'};
  my $gbin = Bio::SeqIO->new(-file=>$$settings{'gb'},-format=>'genbank');
  my $fasta;
  if(defined($$settings{'fasta'})){
    $fasta = Bio::SeqIO->new(-file=>">".$$settings{'fasta'},-format=>'');
    printf STDERR "writing %s... ",$$settings{'fasta'};
  }
  my $gff;
  while( my $seq=$gbin->next_seq()){ #each contig
    my @feats = $seq->all_SeqFeatures(); # features here are loci in the genome
    foreach my $ftr ($seq->all_SeqFeatures()){
      if($ftr->primary_tag eq 'gene'){
        my $seq2fasta = $ftr->seq;
        my $seqid = scalar $seq2fasta->display_name;
        my $locus_tag = ($ftr->get_tag_values('locus_tag'))[0];
        my $seqid_sanitized=$seq->id;
        $seqid_sanitized=~s/_/-/g;
        my $alt_locus_tag=join("_",$seqid_sanitized,$ftr->start,$ftr->end);
        my $defline = sprintf("lcl|%s|%s",$locus_tag,$alt_locus_tag);
        $seq2fasta->display_name($defline); # >lcl|StrainName_Locus|start|end
        my @desc = ();
        if($ftr->has_tag('gene')){push(@desc,$ftr->get_tag_values('gene'));}
        #locate a product tag
        for my $eachftr (@feats){
          if(!$eachftr->has_tag('locus_tag')){next;}
          if($locus_tag ne ($eachftr->get_tag_values('locus_tag'))[0]){next;}
          if($eachftr->primary_tag =~ /CDS|tRNA/){
            if($eachftr->has_tag('product')){
              push(@desc,$eachftr->get_tag_values('product'));
              last;
            }
          }
        }
        if(1>scalar @desc){push(@desc,"predicted cds");}
        $seq2fasta->desc(join(" ",@desc));
        #print STDERR join(" ",$defline,@desc,"\n");
        if($fasta){$fasta->write_seq($seq2fasta) or die "$!\n";}
      }
    }
    #print STDERR "Finished with contig ".$seq->id."\n";
  }
  print STDERR "\ndone\n";
  return 0;
}#end main

sub usage{
  "Usage: ". basename($0) . " --gb=genbank-input-file --fasta=fasta-output-file [--faa=translatedCDS.faa]\n";
}
