#!/usr/bin/env perl

use strict;
use FindBin;
use lib "$FindBin::RealBin/../lib";
$ENV{PATH} = "$FindBin::RealBin:".$ENV{PATH};
use AKUtils qw(logmsg);

use GTTmhmm;

# If TMHMM isn't set, it's the directory one level up from the bin directory
if(!$ENV{TMHMMDIR}){
  $ENV{TMHMMDIR} = `dirname \`which tmhmm\``;
  $ENV{TMHMMDIR}=~s/\/\w+$//i;
  logmsg "Warning: ENV variable for TMHMMDIR is not set.  Setting it here as $ENV{TMHMMDIR}. Avoid this message in the future by adding the following line to ~/.bashrc\n  export TMHMMDIR=$ENV{TMHMMDIR}";
} else {
  logmsg "TMHMMDIR is set as $ENV{TMHMMDIR}";
}

die("Usage: $0 input.faa") unless -f $ARGV[0];

my $inseq = Bio::SeqIO->new(-file => "<$ARGV[0]", -format => 'fasta');

open(TMHMM_OUT_SQL, '>', "$ARGV[0].tmhmm_location.sql") or die;
open(TMHMM_STATS_SQL, '>', "$ARGV[0].tmhmm.sql") or die;

#my $factory = Bio::Tools::Run::Tmhmm->new();
# The TMHMM runner and parser had to be modified to allow parsing of the "inside"/"outside" markers and stats.
my $factory = GTTmhmm->new();
logmsg "Running TMHMM on proteins in $ARGV[0]...";
my $i;
while (my $seq = $inseq->next_seq) {
	$i++; logmsg "Processed $i proteins" if $i % 100 == 0;

	my ($feats, $stats) = $factory->run($seq);
	next if @$feats < 2; # a transmembrane motif results in at least 3 features (inside-transmembrane-outside)

	foreach my $feat (@$feats) {
		my $type = $feat->primary_tag;
		$type =~ s/TMhelix/Transmembrane Helix/;
		print TMHMM_OUT_SQL join("|", ($seq->id, $type, $feat->start, $feat->end))."\n";
	}

	my $tmhmm_measures = ["Length: ", "Number of predicted TMHs: ", "Exp number of AAs in TMHs: ", "Exp number, first 60 AAs: ", "Total prob of N-in: "];
	my @l = ( $seq->id );
	
	for (@$tmhmm_measures) { push(@l, $$stats{$_}); }
	print TMHMM_STATS_SQL join("|", @l)."\n";
}
logmsg "Processed $i proteins, done";

close TMHMM_OUT_SQL;
close TMHMM_STATS_SQL;
