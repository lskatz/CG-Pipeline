#!/usr/bin/env perl
use warnings;
use strict;
use File::Path;
use FindBin;
use lib "$FindBin::RealBin/../lib";
use CGPBase;
my ($organism,$gfffile,$outdir)=@ARGV;
my $fastafile="$outdir/$organism.cds.fna";
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
print "Running formatdb on $fastafile...\n";
makeblastdb($organism,$fastafile,$outdir);
print "BLAST database in $outdir\n";

