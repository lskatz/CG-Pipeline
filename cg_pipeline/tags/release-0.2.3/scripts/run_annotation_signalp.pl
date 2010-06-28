#!/usr/bin/env perl

use strict;
use FindBin;
use lib "$FindBin::RealBin/../lib";
$ENV{PATH} = "$FindBin::RealBin:".$ENV{PATH};
use AKUtils qw(logmsg);

die("Usage: $0 input.faa") unless -f $ARGV[0];

my $seqs = AKUtils::readMfa($ARGV[0]);
my $predictions = AKUtils::getSignalPPredictions($seqs, {signalp_org_type => 'gram-'});

open(SIGNALP_NN_OUT, '>', "$ARGV[0].signalp_nn.sql") or die;
open(SIGNALP_HMM_OUT, '>', "$ARGV[0].signalp_hmm.sql") or die;

foreach my $name (keys %$predictions) {
	foreach my $pred (@{$$predictions{$name}}) {
		if ($$pred{type} eq 'hmm') {
			print SIGNALP_HMM_OUT join("|", ($name, $$pred{decision}, $$pred{probability}, $$pred{CLP}, $$pred{start}, $$pred{end}))."\n";
		} elsif ($$pred{type} eq 'nn') {
			print SIGNALP_NN_OUT join("|", ($name, $$pred{measure}, $$pred{start}, $$pred{end}, $$pred{value}, $$pred{cutoff}, $$pred{decision}))."\n";
		} else { die("Internal error"); }
	}
}

close SIGNALP_NN_OUT;
close SIGNALP_HMM_OUT;
