#!/usr/bin/env perl

my $settings = {};
my $stats = {};

use strict;
use FindBin;
use lib "$FindBin::RealBin/../lib";
use AKUtils qw(logmsg);
use Getopt::Long;
use File::Temp ('tempdir');
use File::Path;
use File::Spec;
use File::Copy;
use File::Basename;
use Bio::SeqIO;
use XML::LibXML::Reader;

$0 = fileparse($0);

my $usage = "$0 [-outfile=parsed_iprscan_output.sql] iprscan_output.xml";

my @cmd_options = ('ipr_hits_out_sql_file=s', 'ipr_classifications_out_sql_file=s', 'ipr_matches_out_sql_file=s', 'ipr_childrefs_out_sql_file=s',);
GetOptions($settings, @cmd_options) or die "Usage: $usage";
die("Usage: $usage") if @ARGV != 1;
#die("Argument blast_db must be supplied") unless $$settings{blast_db};
$$settings{query_xml} = $ARGV[0];
die("File $$settings{query_xml} does not exist") unless -f $$settings{query_xml};

$$settings{ipr_hits_out_sql_file} ||= "$$settings{query_xml}.ipr_hits.sql";
$$settings{ipr_classifications_out_sql_file} ||= "$$settings{query_xml}.ipr_classifications.sql";
$$settings{ipr_matches_out_sql_file} ||= "$$settings{query_xml}.ipr_matches.sql";
$$settings{ipr_childrefs_out_sql_file} ||= "$$settings{query_xml}.ipr_childrefs.sql";

my $file_size = (stat $$settings{query_xml})[7];
logmsg "Processing XML in file $$settings{query_xml}";

my $reader = new XML::LibXML::Reader(location => $$settings{query_xml}) or die "cannot read $$settings{query_xml}\n";
my $parser = XML::LibXML->new();

my (@interpro_hits, @classifications, @matches, @interpro_child_references);
open(IPR_HITS_OUT_SQL, '>', $$settings{ipr_hits_out_sql_file})
	or die "Unable to open file $$settings{ipr_hits_out_sql_file} for writing: $!";
open(IPR_CLASSIFICATIONS_OUT_SQL, '>', $$settings{ipr_classifications_out_sql_file})
	or die "Unable to open file $$settings{ipr_classifications_out_sql_file} for writing: $!";
open(IPR_MATCHES_OUT_SQL, '>', $$settings{ipr_matches_out_sql_file})
	or die "Unable to open file $$settings{ipr_matches_out_sql_file} for writing: $!";
open(IPR_CHILDREFS_OUT_SQL, '>', $$settings{ipr_childrefs_out_sql_file})
	or die "Unable to open file $$settings{ipr_classifications_out_sql_file} for writing: $!";

my $i;
while ($reader->read) {
	if ($reader->name eq 'protein' and $reader->nodeType != XML_READER_TYPE_END_ELEMENT) {
		$i++; logmsg("[".int(100*$reader->byteConsumed/$file_size)."%] Processed $i records...") if $i % 1000 == 0;
		my $entry = $parser->parse_string($reader->readOuterXml);

		my %info;
		$info{protein_id} = $entry->getElementsByTagName('protein')->[0]->attributes->getNamedItem('id')->nodeValue;
		$info{protein_length} = $entry->getElementsByTagName('protein')->[0]->attributes->getNamedItem('length')->nodeValue;

		my @interpro_hits;
		foreach my $hit ($entry->getElementsByTagName('interpro')) {
			my %hit_info;
			$hit_info{query_protein_id} = $info{protein_id};
			$hit_info{query_protein_length} = $info{protein_length};
			$hit_info{interpro_id} = $hit->attributes->getNamedItem('id')->nodeValue;
			$hit_info{interpro_name} = $hit->attributes->getNamedItem('name')->nodeValue;
			$hit_info{interpro_type} = $hit->attributes->getNamedItem('type')->nodeValue;
			push(@interpro_hits, \%hit_info);

			my @l; push(@l, $hit_info{$_}) for qw(query_protein_id query_protein_length interpro_id interpro_name interpro_type);
			s/\|/\\|/g for @l; # escape the pipe characters
			print IPR_HITS_OUT_SQL join('|', @l)."\n";
			$$stats{ipr_hits_record_count}++;

			foreach my $match ($hit->findnodes('match')) {
				my %match_info;
				$match_info{id} = $match->attributes->getNamedItem('id')->nodeValue;
				$match_info{name} = $match->attributes->getNamedItem('name')->nodeValue;
				$match_info{dbname} = $match->attributes->getNamedItem('dbname')->nodeValue;
				foreach my $location ($match->findnodes('location')) {
					my %loc_info;
					foreach my $attr (qw(start end score status evidence)) {
						$loc_info{$attr} = $location->attributes->getNamedItem($attr)->nodeValue;
					}
					$loc_info{protein_id} = $info{protein_id};
					$loc_info{match_id} = $match_info{id};
					$loc_info{match_name} = $match_info{name};
					$loc_info{match_db_name} = $match_info{dbname};
					push(@matches, \%loc_info);

					my @l; push(@l, $loc_info{$_}) for qw(protein_id match_id match_name match_db_name start end score status evidence);
					s/\|/\\|/g for @l; # escape the pipe characters
					print IPR_MATCHES_OUT_SQL join('|', @l)."\n";
					$$stats{ipr_matches_record_count}++;
				}
			}
			next unless $hit->findnodes('child_list');
			foreach my $child_rel_ref ($hit->findnodes('child_list')->[0]->findnodes('rel_ref')) {
				my $ipr_ref = $child_rel_ref->attributes->getNamedItem('ipr_ref')->nodeValue;
				my %child_ref_info = {protein_id => $info{protein_id},
									  interpro_hit_id => $hit_info{interpro_id},
									  interpro_child_reference => $ipr_ref};
				push(@interpro_child_references, \%child_ref_info);

				my @l; push(@l, $child_ref_info{$_}) for qw(protein_id interpro_hit_id interpro_child_reference);
				s/\|/\\|/g for @l; # escape the pipe characters
				print IPR_CHILDREFS_OUT_SQL join('|', @l)."\n";
				$$stats{ipr_childrefs_record_count}++;
			}
		}

		foreach my $c ($entry->getElementsByTagName('classification')) {
			my %classification;
			$classification{query_protein_id} = $info{protein_id};
			$classification{id} = $c->attributes->getNamedItem('id')->nodeValue;
			$classification{class_type} = $c->attributes->getNamedItem('class_type')->nodeValue;
			$classification{category} = $c->findnodes('category')->[0]->firstChild->data;
			$classification{description} = $c->findnodes('description')->[0]->firstChild->data;
			push(@classifications, \%classification);

			my @l; push(@l, $classification{$_}) for qw(query_protein_id id class_type category description);
			s/\|/\\|/g for @l; # escape the pipe characters
			print IPR_CLASSIFICATIONS_OUT_SQL join('|', @l)."\n";
			$$stats{ipr_classifications_record_count}++;

			my @l = $c->findnodes('category'); die if scalar(@l) != 1; @l = $c->findnodes('description'); die if scalar(@l) != 1;
		}
		$reader->next; # skip subtree
	}
}

logmsg "Processed $i records, done";

close IPR_HITS_OUT_SQL;
close IPR_CLASSIFICATIONS_OUT_SQL;
close IPR_MATCHES_OUT_SQL;
close IPR_CHILDREFS_OUT_SQL;

logmsg "Printed $$stats{ipr_hits_record_count} records to $$settings{ipr_hits_out_sql_file}";
logmsg "Printed $$stats{ipr_classifications_record_count} records to $$settings{ipr_classifications_out_sql_file}";
logmsg "Printed $$stats{ipr_matches_record_count} records to $$settings{ipr_matches_out_sql_file}";
logmsg "Printed $$stats{ipr_childrefs_record_count} records to $$settings{ipr_childrefs_out_sql_file}";

#my $io = Bio::SeqIO->new(-format => "interpro", -file => $ARGV[0]);
#while (my $seq = $io->next_seq) { print "seq:\n\tid:".$seq->accession."\n"; }
