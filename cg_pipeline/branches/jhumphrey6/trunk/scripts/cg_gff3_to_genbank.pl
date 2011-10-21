#!/usr/bin/env perl
# author: Jay Humphrey <jhumphrey6@gmail.com>
# convert gff to genbank

use warnings;
use strict;
use File::Path;
use FindBin;
use lib "$FindBin::RealBin/../lib";
#use CGPipeline::Store; 
use CGPipeline::Store::Prediction; 
$0 =~ s/^.*\/([^\/]+)$/$1/;
if(!@ARGV){print "usage: $0 input_fasta_file input_gff_file output_genbank_filename\n";exit;}
my ($fnain,$gffin,$gbout)=@ARGV;
#my $sqlitedb=sprintf("%s.%d.sqlite",$gbout,time);
my $sqlitedb=sprintf("%s.sqlite",$gbout);
print "Creating $sqlitedb...\n";
my $cgp = CGPipeline::Store::Prediction->new(file=>$sqlitedb,verbose=>1,volatile=>1,locus_id_format=>'M16207_loc_%04d'); printf "object %s\n",ref $cgp;
$cgp->load_fasta($fnain);
$cgp->clean;
$cgp->load_gff($gffin);
$cgp->write_genbank($gbout);
$cgp = undef;
print "Removing $sqlitedb...\n";
#system("rm $sqlitedb");
