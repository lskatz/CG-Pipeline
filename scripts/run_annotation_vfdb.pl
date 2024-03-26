#!/usr/bin/env perl

# run-assembly: Perform standard assembly protocol operations on 454 pyrosequenced flowgram file(s)
# Author: Andrey Kislyuk (kislyuk@gatech.edu)

package PipelineRunner;
my ($VERSION) = ('$Id: $' =~ /,v\s+(\d+\S+)/o);

my $settings = {
	appname => 'cgpipeline',
	report_top_vfdb_hits => 1,
};
my $stats;

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
use List::Util qw(min max sum shuffle);
use CGPipelineUtils;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
use Data::Dumper;
$|++;

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };

exit(main());

sub main() {
	$settings = AKUtils::loadConfig($settings);

	my @cmd_options = qw(outfile=s blastfile=s db=s);
	GetOptions($settings, @cmd_options) or die;

  $$settings{vfdb_blast_db}||=File::Spec->rel2abs($$settings{db});

	die("Usage: $0 input.mfa") if @ARGV != 1;
	for (qw(vfdb_blast_db min_vfdb_aa_coverage min_vfdb_aa_identity)) {
		die("Argument $_ must be supplied") unless $$settings{$_};
	}
	$$settings{query_mfa} = $ARGV[0];
	die("File $$settings{query_mfa} does not exist") unless -f $$settings{query_mfa};
	$$settings{outfile} ||= "$$settings{query_mfa}.vfdb_hits.sql";

	$$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($$settings{keep}));
	logmsg "Temporary directory is $$settings{tempdir}";

  # don't blast if a blastfile has been given
  my $report;
  if(!$$settings{blastfile}){
    $ENV{BLASTDB} = (fileparse($$settings{vfdb_blast_db}))[1];
    my $bf = Bio::Tools::Run::StandAloneBlast->new(database => $$settings{vfdb_blast_db},
      outfile => "$$settings{tempdir}/$0.$$.blast_out",
      program => 'blastp',
      a => AKUtils::getNumCPUs());

    logmsg "Running BLAST on $$settings{query_mfa} vs. $$settings{vfdb_blast_db}...";
    $report = $bf->blastall($$settings{query_mfa});
  } else {
    $report=new Bio::SearchIO(-format=>'blastxml',-file=>$$settings{blastfile});
  }

	my $query_seqs = AKUtils::readMfa($$settings{query_mfa});
#	my $db_seqs = AKUtils::readMfa($$settings{db_mfa});

	my %reported_hits;

	open(OUT, '>', $$settings{outfile}) or die;
	while (my $result = $report->next_result) {
		while (my $hit = $result->next_hit) {
			HSP: while (my $hsp = $hit->next_hsp) {
				my ($id1, $id2) = ($hsp->seq('hit')->id, $hsp->seq('query')->id);
				die("Internal error - $id1 hit against itself") if $id1 eq $id2;
				my $l2 = length($$query_seqs{$id2}) || $hsp->seq->length; die("Internal error - something wrong with the query seq") unless $l2;
				my $q_coverage = $hsp->length('query') / $l2;
#				my $t_coverage = $hsp->length('hit') / $l1;
				die("Internal error - coverage is > 100%") if $q_coverage > 1; # or $t_coverage > 1;

				if ($q_coverage > $$settings{min_vfdb_aa_coverage}
#					and $t_coverage > $$settings{min_vfdb_aa_coverage}
          and $hsp->frac_conserved >= $$settings{min_vfdb_aa_similarity}
					and $hsp->percent_identity >= $$settings{min_vfdb_aa_identity} * 100) {
                    $reported_hits{$id2}->{$id1} = {query_id => $id2, target_id => $id1, evalue => $hsp->evalue,
													query_coverage => sprintf("%.2f", $q_coverage), db_name => $$settings{vfdb_blast_db},
													percent_identity => sprintf("%.2f", $hsp->percent_identity)};
				}
			}
		}
	}

	foreach my $query_id (keys %reported_hits) {
		my @sorted_hits = sort {$$a{evalue} <=> $$b{evalue}} values(%{$reported_hits{$query_id}});
		foreach my $hit (@sorted_hits[0..$$settings{report_top_vfdb_hits}-1]) {
			my @l;
			push(@l, $$hit{$_}) for qw(query_id target_id evalue query_coverage db_name percent_identity);
			s/\|/\\|/g for @l; # escape the pipe characters
			print OUT join('|', @l)."\n";
		}
	}
			
	close OUT;
	logmsg "Report is in $$settings{outfile}";
	return 0;
}
