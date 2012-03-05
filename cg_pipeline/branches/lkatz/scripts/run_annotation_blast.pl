#!/usr/bin/env perl

# run-assembly: Perform standard assembly protocol operations on 454 pyrosequenced flowgram file(s)
# Author: Andrey Kislyuk (kislyuk@gatech.edu)

package PipelineRunner;
my ($VERSION) = ('$Id: $' =~ /,v\s+(\d+\S+)/o);

my $settings = {
	appname => 'cgpipeline',
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

	my @cmd_options = qw(outfile=s blastfile=s db=s parametersForBlast=s min_aa_coverage=i min_aa_identity=i min_aa_similarity=i tempdir=s keep);
	GetOptions($settings, @cmd_options) or die;

  $$settings{blast_db}||=File::Spec->rel2abs($$settings{db});

	die(usage()) if @ARGV != 1;
  $$settings{min_aa_coverage}||=10;
  $$settings{min_aa_identity}||=10;
  $$settings{min_aa_similarity}||=10;
	for (qw(blast_db)) {
		die("Argument $_ must be supplied") unless $$settings{$_};
	}
	$$settings{query_mfa} = $ARGV[0];
	die("File $$settings{query_mfa} does not exist") unless -f $$settings{query_mfa};
	$$settings{outfile} ||= "$$settings{query_mfa}.unspecifiedhits.sql";

	$$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($$settings{keep}));
	logmsg "Temporary directory is $$settings{tempdir}";

  # don't blast if a blastfile has been given
  my $report;
  if(!$$settings{blastfile}){
    $ENV{BLASTDB} = (fileparse($$settings{blast_db}))[1];
    my $numcpus=AKUtils::getNumCPUs();
    my $command="legacy_blast.pl blastall -p blastp -a $numcpus -o $$settings{tempdir}/$0.$$.blast_out -d $$settings{blast_db} -i $$settings{query_mfa} $$settings{parametersForBlast} -m 7";
    logmsg "Running BLAST on $$settings{query_mfa} vs. $$settings{blast_db}...\n  $command";
    system($command);
    die "Problem with blast" if $?;
    logmsg "Finished with BLAST";
    $report=new Bio::SearchIO(-format=>'blastxml',-file=>"$$settings{tempdir}/$0.$$.blast_out");
  } else {
    $report=new Bio::SearchIO(-format=>'blastxml',-file=>$$settings{blastfile});
  }

  my $resultCounter=0;
  logmsg "Reading $$settings{query_mfa} and parsing blast result";
	my $query_seqs = AKUtils::readMfa($$settings{query_mfa});

	my %reported_hits;

	open(OUT, '>', $$settings{outfile}) or die;
	while (my $result = $report->next_result) {
		while (my $hit = $result->next_hit) {
			HSP: while (my $hsp = $hit->next_hsp) {
				my ($id1, $id2) = ($hsp->seq('hit')->id, $hsp->seq('query')->id);
				die("Internal error - $id1 hit against itself") if $id1 eq $id2;
				my $l2 = length($$query_seqs{$id2}) || $hsp->seq->length; die("Internal error - something wrong with the query seq") unless $l2;
				my $q_coverage = $hsp->length('query') / $l2 * 100;
#				my $t_coverage = $hsp->length('hit') / $l1;
        my $percent_conserved=$hsp->frac_conserved*100;
        my $percent_identity=$hsp->percent_identity;
        #print join("\t",$id1,$id2,"$q_coverage>$$settings{min_aa_coverage}","$percent_conserved>= $$settings{min_aa_similarity}","$percent_identity>= $$settings{min_aa_identity}")."\n"; # debug
				die("Internal error - coverage is > 100%") if $q_coverage > 100; 

				if ($q_coverage > $$settings{min_aa_coverage}
          and $percent_conserved >= $$settings{min_aa_similarity}
					and $percent_identity  >= $$settings{min_aa_identity}) {
              my (undef, $hit_accession, $hit_name) = split(/\|/, $hit->name);
              $reported_hits{$id2}->{$id1} = {
                query_id => $id2, target_id => $id1, evalue => $hsp->evalue,
                query_coverage => sprintf("%.2f", $q_coverage), db_name => $$settings{blast_db},
                percent_identity => sprintf("%.2f", $percent_identity),
                hspLength=>$hit->length,
                description=>$hit->description || '.',
                rank=>$hit->rank,
                score=>$hit->score,
                bits=>$hit->bits,
                percent_conserved=>$percent_conserved,
                hit_accession=>$hit_accession,
                hit_name=>$hit_name,
              };
				}
			}
		}
    $resultCounter++; logmsg "Processed $resultCounter proteins" if($resultCounter % 100 == 0);
	}

	foreach my $query_id (keys %reported_hits) {
		my @sorted_hits = sort {$$a{evalue} <=> $$b{evalue}} values(%{$reported_hits{$query_id}});
		foreach my $hit (@sorted_hits) {
			my @l;
			push(@l, $$hit{$_}) for qw(query_id target_id evalue query_coverage db_name percent_identity hspLength description rank score bits percent_conserved hit_accession hit_name);
			s/\|/\\|/g for @l; # escape the pipe characters
			print OUT join('|', @l)."\n";
		}
	}
			
	close OUT;
	logmsg "Report is in $$settings{outfile}";
	return 0;
}

sub usage{
  "Usage: $0 input.mfa
  -d blastdatabase
    blast formatted database
  -p 'blast parameters'
    A quoted string to pass to blastall (legacy_blast.pl blastall -p blastp is used internally)
    There is no sanity check for this parameter.
  -o outfile.sql
    pipe delimited output file

  -b blast results file
    Optionally to skip blast and use your own blast file
  --min_aa_coverage integer 1 to 100
  --min_aa_identity integer 1 to 100
  --min_aa_similarity integer 1 to 100
  -t tempdir (default is a subdir in /tmp)
  "
}
