#!/usr/bin/env perl

# run-annotation-isElements: annotate IS elements from aa.fasta file
# Author: Lee Katz <lkatz@cdc.gov>

package PipelineRunner;
my ($VERSION) = ('$Id: $' =~ /,v\s+(\d+\S+)/o);

my $settings = {
	appname => 'cgpipeline',
  is_blast_db=>"$ENV{HOME}/bin/cg_pipeline/data/is.faa",
  is_min_aa_coverage=>0.5,
  is_min_aa_identity=>0.3,
  is_max_aa_evalue=>1E-5,
  is_report_top_hits=>10,
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

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };

exit(main());

sub main() {
	$settings = AKUtils::loadConfig($settings);
	die(usage()) if @ARGV < 1;

	my @cmd_options = ('outfile=s','tempdir=s','infile=s');
	GetOptions($settings, @cmd_options) or die;

	for (qw(is_blast_db is_min_aa_coverage is_min_aa_identity is_max_aa_evalue)) {
		die("Argument $_ must be supplied") unless $$settings{$_};
	}
	$$settings{query_mfa} ||= $$settings{infile} or die usage();
	die("File $$settings{query_mfa} does not exist") unless -f $$settings{query_mfa};
	$$settings{outfile} ||= "$$settings{query_mfa}.is_hits.sql";

	$$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($$settings{keep}));
	logmsg "Temporary directory is $$settings{tempdir}";

  whichBlast($settings); # which blast? blast+?  Modify the settings.
  $$settings{num_cpus}||=AKUtils::getNumCPUs();

	$ENV{BLASTDB} = (fileparse($$settings{is_blast_db}))[1];
	my $bf = Bio::Tools::Run::StandAloneBlast->new(database => $$settings{is_blast_db},
		outfile => "$$settings{tempdir}/$0.$$.blast_out",
		program => 'blastp',
		a => $$settings{num_cpus},
    e => $$settings{is_max_aa_evalue},
    F => 'F', # filter off
    W => '2', # reduced word size
  );

	logmsg "Running BLAST on $$settings{query_mfa} vs. $$settings{is_blast_db}...";
  # blast+ query
  my $report;
  if($$settings{blast_version} eq 'blast+'){
    die "BLAST+ is not supported right now";
    # filtering for blastp is off by default
    my $blastCommand="blastp -query $$settings{query_mfa} -db $$settings{is_blast_db} -num_threads $$settings{num_cpus} -out $$settings{tempdir}/$0.$$.blast_out -evalue=>$$settings{is_max_aa_evalue} -word_size 2 ";
    system($blastCommand); die if $?;
    $report=new Bio::SearchIO(-format=>'blast',-file=>"$$settings{tempdir}/$0.$$.blast_out");
  }
  # regular blast query
  else {
    $report=$bf->blastall($$settings{query_mfa});
  }

	my $query_seqs = AKUtils::readMfa($$settings{query_mfa});

	my %reported_hits;

	open(OUT, '>', $$settings{outfile}) or die;
	while (my $result = $report->next_result) {
		while (my $hit = $result->next_hit) {
			HSP: while (my $hsp = $hit->next_hsp) {
				my ($id1, $id2) = ($hsp->seq('hit')->id, $hsp->seq('query')->id);
				die if $id1 eq $id2;
				my $l2 = length($$query_seqs{$id2}); die("Internal error") unless $l2;
				my $q_coverage = $hsp->length('query') / $l2;
				die("Internal error") if $q_coverage > 1; # or $t_coverage > 1;

				if ($q_coverage > $$settings{is_min_aa_coverage}
					and $hsp->percent_identity >= $$settings{is_min_aa_identity} * 100) {
                    $reported_hits{$id2}->{$id1} = {query_id => $id2, target_id => $id1, evalue => $hsp->evalue,
													query_coverage => sprintf("%.2f", $q_coverage), db_name => $$settings{is_blast_db},
													percent_identity => sprintf("%.2f", $hsp->percent_identity)};
				}
			}
		}
	}

	foreach my $query_id (keys %reported_hits) {
		my @sorted_hits = sort {$$a{evalue} <=> $$b{evalue}} values(%{$reported_hits{$query_id}});
    
    splice(@sorted_hits,$$settings{is_reported_top_hits}+1) if(@sorted_hits>$$settings{is_reported_top_hits});
		foreach my $hit (@sorted_hits) {
			my @l;
			push(@l, $$hit{$_}) for qw(query_id target_id evalue query_coverage db_name percent_identity);
			s/\|/\\|/g for @l; # escape the pipe characters
			print OUT join('|', @l)."\n";
		}
	}
			
	close OUT;

  # end stats
  my $numIsElements=`wc -l $$settings{outfile} | cut -f 1 | sort|uniq` + 0;
	logmsg "$numIsElements IS elements found. Report is in $$settings{outfile}";
	return 0;
}

# modify settings to figure out which version of blast we are dealing with
# TODO if need be, put this into AKUtils
# TODO allow for whichBlast to pick up on the other four flavors of blast
sub whichBlast{
  my($settings)=@_;
  my $version="before2009";
  $$settings{path_to_blast}=AKUtils::fullPathToExec("blastall")||AKUtils::fullPathToExec("blastp")||"";
  if($$settings{path_to_blast}=~/blastall/){
    $$settings{blast_version}="before2009";
  }
  elsif($$settings{path_to_blast}=~/blastp/){
    $$settings{blast_version}="blast+";
  }
}

sub usage{
"
This script finds IS elements in a protein multi fasta file
Usage: $0 -i input.mfa
  -i input multifasta file (protein only)
  -o output file
    default: input.mfa.is_hits.sql (depending on what your infile is named)
  -t temporary directory
"
}
