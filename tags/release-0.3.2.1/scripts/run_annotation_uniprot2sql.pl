#!/usr/bin/env perl

# Take the output of run_annotation_blast.pl and print data for Uniprot IDs mentioned, using the database made by parse_uniprot_xml.pl
# Author: Andrey Kislyuk

my $settings = {};

use strict;
use FindBin;
use lib "$FindBin::RealBin/../lib";
$ENV{PATH} = "$FindBin::RealBin:".$ENV{PATH};
use AKUtils qw(logmsg);

use Getopt::Long;
use File::Path;
use File::Spec;
use File::Copy;
use File::Basename;
use BerkeleyDB;

my $usage = "$0 -uniprot_db3=u.db3 -uniprot_evidence_db3=ue.db3 [-uniprot_sql_outfile=u.sql] [-uniprot_evidence_sql_outfile=e.sql] blast.sql";

my @cmd_options = ('blast_db=s', 'uniprot_db3=s', 'uniprot_evidence_db3=s', 'uniprot_sql_outfile=s', 'uniprot_evidence_sql_outfile=s');

GetOptions($settings, @cmd_options) or die "Usage: $usage";
die("Usage: $usage") if @ARGV != 1;
die("Argument uniprot_db3 must be supplied") unless $$settings{uniprot_db3};
die("Argument uniprot_evidence_db3 must be supplied") unless $$settings{uniprot_evidence_db3};
$$settings{blast_sql_file} = $ARGV[0];
die("File $$settings{blast_sql_file} does not exist") unless -f $$settings{blast_sql_file};

my (%uniprot_h, %uniprot_evidence_h);
tie(%uniprot_h, "BerkeleyDB::Hash", -Filename => $$settings{uniprot_db3}, -Flags => DB_RDONLY, -Property => DB_DUP)
	or die "Cannot open file $$settings{uniprot_db3}: $! $BerkeleyDB::Error\n";
tie(%uniprot_evidence_h, "BerkeleyDB::Hash", -Filename => $$settings{uniprot_evidence_db3}, -Flags => DB_RDONLY, -Property => DB_DUP)
	or die "Cannot open file $$settings{uniprot_evidence_db3}: $! $BerkeleyDB::Error\n";

$$settings{uniprot_sql_outfile} ||= "$$settings{blast_sql_file}.uniprot.sql";
$$settings{uniprot_evidence_sql_outfile} ||= "$$settings{blast_sql_file}.uniprot_evidence.sql";

open(UNIPROT_OUT_SQL, '>', $$settings{uniprot_sql_outfile})
	or die "Unable to open file $$settings{uniprot_sql_outfile} for writing: $!";
open(UNIPROT_EVIDENCE_OUT_SQL, '>', $$settings{uniprot_evidence_sql_outfile})
	or die "Unable to open file $$settings{uniprot_evidence_sql_outfile} for writing: $!";

logmsg "Extracting Uniprot records for BLAST results in $$settings{blast_sql_file}";

my $i;
open(IN, '<', $$settings{blast_sql_file}) or die;
while (<IN>) {
	chomp;
	my (undef, $accession) = split /\|/;
	my $uniprot_line = $uniprot_h{$accession};
	my $uniprot_evidence_line = $uniprot_evidence_h{$accession};
	die("No record found for accession $accession") unless $uniprot_line;

	print UNIPROT_OUT_SQL $uniprot_line."\n";
	print UNIPROT_EVIDENCE_OUT_SQL $uniprot_evidence_line."\n";
	$i++; logmsg "Processed $i records..." if $i % 1000 == 0;
}
logmsg "Processed $i records, done";
close IN;

close UNIPROT_OUT_SQL;
close UNIPROT_EVIDENCE_OUT_SQL;
