#!/usr/bin/env perl
# author: Jay Humphrey <jhumphrey6@gatech.edu>
# Creates a database from a given annotated gff file

use warnings;
use strict;
use File::Path;
use FindBin;
use lib "$FindBin::RealBin/../lib";
use CGPBase;
if(!@ARGV){print "$0 organism output-gff-file output-directory\n";exit;}
my ($organism,$gfffile,$outdir)=@ARGV;
my $fastafile="$outdir/$organism.cds.fna";
my $fastaprotfile="$outdir/$organism.cds.faa";
my $sqlitedb="$outdir/$organism.sqlite";
if(! -d $outdir ){
	print "Creating folder $outdir...\n";	
	mkpath([$outdir]) or die "$!\n";
}
if( -e $sqlitedb){
	print "Removing old $sqlitedb...\n";	
	system("rm $sqlitedb");
}
print "Creating $sqlitedb...\n";
create_db($sqlitedb,$gfffile);
if( -e $fastafile){
	print "Removing old $fastafile...\n";	
	system("rm $fastafile");
}
print "Creating $fastafile...\n";
sqlite2fasta($sqlitedb,$fastafile);
sqlite2fastaprot($sqlitedb,$fastaprotfile);
print "Running formatdb on $fastafile...\n";
makeblastdb($organism,$fastafile,'F',$outdir);
makeblastdb($organism,$fastaprotfile,'T',$outdir);
print "BLAST database in $outdir\n";

