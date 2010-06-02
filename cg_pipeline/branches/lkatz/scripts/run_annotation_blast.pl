#!/usr/bin/env perl

my $settings = {};
my $stats = {};

use strict;
use FindBin;
use lib "$FindBin::RealBin/../lib";
$ENV{PATH} = "$FindBin::RealBin:".$ENV{PATH};
use AKUtils qw(logmsg);

use Getopt::Long;
use File::Temp ('tempdir');
use File::Path;
use File::Spec;
use File::Copy;
use File::Basename;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;

$0 = fileparse($0);

my $usage = "$0 -blast_db=db [-outfile=blast.sql] input.mfa";

my @cmd_options = ('blast_db=s', 'tempdir=s', 'outfile=s', 'evidence_outfile=s', 'keep');
GetOptions($settings, @cmd_options) or die "Usage: $usage";
die("Usage: $usage") if @ARGV != 1;
die("Argument blast_db must be supplied") unless $$settings{blast_db};
$$settings{query_mfa} = $ARGV[0];
die("File $$settings{query_mfa} does not exist") unless -f $$settings{query_mfa};

$$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($$settings{keep}));
logmsg "Temporary directory is $$settings{tempdir}";

my $bf = Bio::Tools::Run::StandAloneBlast->new(database => $$settings{blast_db},
											   outfile => "$$settings{tempdir}/$0.$$.blast_out",
											   program => 'blastp',
											   a => AKUtils::getNumCPUs());

my $inseq = Bio::SeqIO->new(-file => "<$$settings{query_mfa}", -format => 'fasta');

logmsg "Running BLAST on $$settings{query_mfa} vs. [$$settings{blast_db}]...";
#my $report = $bf->blastall($$settings{query_mfa});

$$settings{blast_sql_file} = $$settings{outfile};
$$settings{blast_sql_file} ||= "$$settings{query_mfa}.blast.sql";
open(BLAST_OUT_SQL, '>', $$settings{blast_sql_file}) or die "Unable to open file $$settings{blast_sql_file} for writing: $!";

my $i;
while (my $seq = $inseq->next_seq) {
	my $report = $bf->blastall($seq);
	my $h;
	while (my $result = $report->next_result) {
		$h++; last if $h > 15;
		while (my $hit = $result->next_hit) {
			next unless $hit->{_hsps}; # bug in Bio::Search::Hit::GenericHit?
			my $hsp = $hit->next_hsp; # TODO: verify that the first hsp returned is always ranked best
			my (undef, $hit_accession, $hit_name) = split(/\|/, $hit->name);

			my @l = ($seq->id, $hit_accession, $hit->length, $hit->description, $hit->rank,
					 $hsp->score, $hsp->bits, $hsp->evalue, int($hsp->frac_identical * 100), int($hsp->frac_conserved * 100));
			s/\|/\\|/g for @l; # escape the pipe characters

			print BLAST_OUT_SQL join('|', @l)."\n";
			$$stats{hit_record_count}++;
		}
	}
	$i++; logmsg "Processed $i proteins" if $i % 100 == 0;
}
logmsg "Processed $i proteins, done";

close BLAST_OUT_SQL;

logmsg "Printed $$stats{hit_record_count} hits to $$settings{blast_sql_file}";
