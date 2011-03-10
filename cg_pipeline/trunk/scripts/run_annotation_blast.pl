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

whichBlast($settings); # which blast? blast+?  Modify the settings.
$$settings{num_cpus}=AKUtils::getNumCPUs();
$$settings{blast_max_hits} ||= 250;

# TODO speed up blast by limiting the number of results
my $bf = Bio::Tools::Run::StandAloneBlast->new(database => $$settings{blast_db},
											   outfile => "$$settings{tempdir}/$0.$$.blast_out",
											   program => 'blastp',
											   a => $$settings{num_cpus});


my $inseq = Bio::SeqIO->new(-file => "<$$settings{query_mfa}", -format => 'fasta');

logmsg "Running BLAST on $$settings{query_mfa} vs. [$$settings{blast_db}]...";
#my $report = $bf->blastall($$settings{query_mfa});

$$settings{blast_sql_file} = $$settings{outfile};
$$settings{blast_sql_file} ||= "$$settings{query_mfa}.blast.sql";
open(BLAST_OUT_SQL, '>', $$settings{blast_sql_file}) or die "Unable to open file $$settings{blast_sql_file} for writing: $!";

my $i;
while (my $seq = $inseq->next_seq) {
  my $report;
  
  # perform the query
  if($$settings{blast_version} eq 'blast+'){
    my $seqString=">".$seq->id."\n".$seq->seq;
    my $blastCommand="echo '$seqString'|blastp -query -db $$settings{blast_db} -num_threads $$settings{num_cpus} -out $$settings{tempdir}/$0.$$.blast_out ";
    #$blastCommand.="-num_descriptions $$settings{blast_max_hits} -num_alignments $$settings{blast_max_hits} ";
    system($blastCommand); die if $?;
    #system("echo $$settings{tempdir}/$0.$$.blast_out");exit;
    $report=new Bio::SearchIO(-format=>'blast',-file=>"$$settings{tempdir}/$0.$$.blast_out");
  }
  else{
  	$report = $bf->blastall($seq);
  }

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

