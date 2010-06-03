#!/usr/bin/env perl

# "Naive" combination and layout of Newbler and AMOScmp mapped assemblies

# Usage: prep-hints.pl 454AlignmentInfo.tsv amos-layout.layout > hints.fasta

use strict;
use AKUtils 'intersectLength';

die if @ARGV < 4;
my ($newbler_aligns, $newbler_seqs, $amos_layout, $amos_seqs) = @ARGV;
my $amos_contig_seqs = AKUtils::readMfa($amos_seqs);
my $newbler_contig_seqs = AKUtils::readMfa($newbler_seqs);

my (%newbler_contigs, %newbler_gaps, %amos_contigs, %amos_gaps);

open(FH, '<', $newbler_aligns) or die;
my ($cur_seq, $cur_start, $cur_end);
while (<FH>) {
	if (/^\>([^\s]+)\s+(\d+)$/) {
		if ($cur_end) {
			push(@{$newbler_contigs{$cur_seq}}, [$cur_start, $cur_end]);
			push(@{$newbler_gaps{$cur_seq}}, {seq => $cur_seq, start => $cur_end+1, end => $2-1, length => $2-$cur_end}) if $cur_seq eq $1;
		}
		($cur_seq, $cur_start) = ($1, $2);
	} else {
		/^(\d+)\s/;
		$cur_end = $1;
	}
}
push(@{$newbler_contigs{$cur_seq}}, [$cur_start, $cur_end]);
close FH;

open(FH, '<', $amos_layout) or die;
my $c;
while (<FH>) {
	if (/^C\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([\d\-]+)-([\d\-]+)$/) {
		$c++;
		die("$1 != $c") if $1 != $c;
		$amos_contigs{$3}->{$4} = {contig_id => $1, seq => $3, start => $4, end => $5, length=>$5-$4+1};
	}
}
close FH;

foreach my $seq (keys %amos_contigs) {
	my $prev_end = 0;
	foreach my $start (sort {$a <=> $b} keys %{$amos_contigs{$seq}}) {
		if ($start > $prev_end) {
			push(@{$amos_gaps{$seq}}, {seq => $seq, start => $prev_end+1, end => $start-1, length => $start-$prev_end-1});
		}
		$prev_end = $amos_contigs{$seq}->{$start}->{end};
	}
}

my %gap_closer_ids;
my @unclosed_gaps;
my $min_overhang = 40;
my $out_contig_id = 1;
my ($unclosed_gap_length, $partial_gap_length);
# if contig spans both overhangs, admit on length
# 
foreach my $seq (keys %newbler_gaps) {
	my $i;
    GAP: foreach my $gap (@{$newbler_gaps{$seq}}) {
		$i++;
		print "$i\t$$gap{start}\t$$gap{length}\n";

		my ($best_contig, $best_ol, $best_oh) = (undef, -1e999, -1e999);
		foreach my $contig (values %{$amos_contigs{$seq}}) {
			if ((AKUtils::intersectLength($$gap{start}, $$gap{end}, $$contig{start}, $$contig{end}) > $min_overhang)
				or (AKUtils::intersectLength($$gap{start}, $$gap{end}, $$contig{start}, $$contig{end}) >= $$gap{end}-$$gap{start})) {
				my ($overhang_lo, $overhang_hi) = ($$gap{start}-$$contig{start}, $$contig{end}-$$gap{end});
				if ($overhang_lo > $min_overhang and $overhang_hi > $min_overhang) {
					$best_contig ||= $contig;
					if ($$best_contig{length} < $$contig{length}) {
						$best_contig = $contig; $best_ol = $overhang_lo; $best_oh = $overhang_hi;
					}
				} elsif ($overhang_lo > $best_ol and $overhang_hi > $best_oh) {
					$best_contig = $contig; $best_ol = $overhang_lo; $best_oh = $overhang_hi;
				}
			}
		}
		if ($best_contig) {
=head1
			my $interval_length = 1000;
			if ($$best_contig{length} < $interval_length) {
				$gap_closer_ids{$$best_contig{contig_id}} = 1;
				next GAP;
			} else { # need to chop it up
				my $start_in_contig = $$gap{start} - ($min_overhang*4) - $$best_contig{start}; $start_in_contig = 1 if $start_in_contig < 1;
				my $length_in_contig = $$gap{length} + ($min_overhang*8);
				
				for (my $pos_in_contig = $start_in_contig;
					 $pos_in_contig < $start_in_contig+$length_in_contig+$interval_length;
					 $pos_in_contig += $interval_length - (2*$min_overhang)) {
					last if $pos_in_contig > $$best_contig{length}; # um...
					print ">bridge_contig".$out_contig_id." $seq:$$best_contig{start}..$$best_contig{end}:$pos_in_contig\n";
					my $adj_pos_in_contig = $pos_in_contig;
					$adj_pos_in_contig = $$best_contig{length} - $interval_length if $adj_pos_in_contig + $interval_length > $$best_contig{length};
					my $seq = substr($$amos_contig_seqs{$$best_contig{contig_id}}, $adj_pos_in_contig, $interval_length);
					$seq =~ s/(.{80})/$1\n/g;
					$seq .= "\n" unless $seq =~ /\n$/;
					print "$seq";
					$out_contig_id++;
				}
			}
=cut
#                print "$$best_contig{start}..$$best_contig{end} contains $$gap{start}..$$gap{end}, overhangs $best_ol, $best_oh, sl = $start_in_contig..$length_in_contig\n";
            $partial_gap_length += ($$gap{end}-$$gap{start}) - intersectLength($$gap{start}, $$gap{end}, $$best_contig{start}, $$best_contig{end});
		} else {
			$unclosed_gap_length += $$gap{length};
			push(@unclosed_gaps, $gap);
		}
	}
}

foreach my $gap (@unclosed_gaps) { print "Unclosed gap: $$gap{start}\t$$gap{length}\n"; }
print "Total unclosed gaps: ".@unclosed_gaps.", length $unclosed_gap_length (partial $partial_gap_length)\n";

__END__
foreach my $id (sort {$a<=>$b} keys %gap_closer_ids) {
	print ">$id\n";
	my $seq = $$amos_contig_seqs{$id};
	$seq =~ s/(.{80})/$1\n/g;
	$seq .= "\n" unless $seq =~ /\n$/;
	print "$seq";
}


#foreach my $seq (keys %amos_gaps) {
#	foreach my $gap (@{$amos_gaps{$seq}}) {
#		print "AMOS: $$gap{start}\t$$gap{end}\t$$gap{length}\n";
#	}
#}
