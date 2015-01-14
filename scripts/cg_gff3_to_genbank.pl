#!/usr/bin/env perl
# author: Jay Humphrey <jhumphrey6@gatech.edu>
# convert gff to genbank

use warnings;
use strict;
use File::Path;
use FindBin;
use lib "$FindBin::RealBin/../lib";
use CGPipeline::SQLiteDB; 
$0 =~ s/^.*\/([^\/]+)$/$1/;
if(!@ARGV){print "usage: $0 input_fasta_file input_gff_file output_genbank_filename\n";exit;}
my ($fnain,$gffin,$gbout)=@ARGV;
my $sqlitedb=sprintf("%s.%d.sqlite",$gbout,time);
print "Creating $sqlitedb...\n";
my $cgp = CGPipeline::SQLiteDB->new(file=>$sqlitedb);
$cgp->load_fasta($fnain);
$cgp->load_gff($gffin);
$cgp->write_genbank($gbout);
print "Removing $sqlitedb...\n";
system("rm $sqlitedb");
