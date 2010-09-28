#!/usr/bin/env perl

# Parse uniprot XML file and output DB3 file with concise information that we need.
# Author: Andrey Kislyuk

use strict;
use FindBin;
use lib "$FindBin::RealBin/../lib";
$ENV{PATH} = "$FindBin::RealBin:".$ENV{PATH};
use AKUtils qw(logmsg);

use XML::LibXML::Reader;
use BerkeleyDB;
#use threads;
#use Thread::Queue;

die("Usage: $0 uniprot_sprot.xml uniprot_trembl.xml") if @ARGV < 1;
my @infiles = @ARGV;
for (@infiles) { die("File $_ not found") unless -f $_; }

my (%uniprot_h, %uniprot_evidence_h);
my ($dbfile, $evidence_dbfile) = ("cgpipeline.db3", "cgpipeline.evidence.db3");
unlink $dbfile; unlink $evidence_dbfile;
tie(%uniprot_h, "BerkeleyDB::Hash", -Filename => $dbfile, -Flags => DB_CREATE, -Property => DB_DUP)
	or die "Cannot open file $dbfile: $! $BerkeleyDB::Error\n";
tie(%uniprot_evidence_h, "BerkeleyDB::Hash", -Filename => $evidence_dbfile, -Flags => DB_CREATE, -Property => DB_DUP)
	or die "Cannot open file $evidence_dbfile: $! $BerkeleyDB::Error\n";

my $j;
foreach my $infile (@infiles) {
	logmsg "Processing XML in file $infile...";

	my $file_size = (stat $infile)[7];

	my $reader = new XML::LibXML::Reader(location => $infile) or die "cannot read $infile\n";	

	my $i;
	while ($reader->read) {
		if ($reader->name eq 'entry' and $reader->nodeType != XML_READER_TYPE_END_ELEMENT) {
			$i++; $j++; logmsg("[".int(100*$reader->byteConsumed/$file_size)."%] Processed $i records...") if $i % 1000 == 0;
			processUniprotXMLEntry($reader->readOuterXml);
			$reader->next; # skip subtree
		}
	}
	logmsg "Processed $i records, done with file $infile";
}
logmsg "Processed $j records, done";

untie %uniprot_h;
untie %uniprot_evidence_h;

sub processUniprotXMLEntry($) {
	my ($xml_string) = @_;
	my $parser = XML::LibXML->new();
	my $entry = $parser->parse_string($xml_string);

	my %info;
	$info{accession} = $entry->getElementsByTagName('accession')->[0]->firstChild->nodeValue;

	$info{dataset} = $entry->getElementsByTagName('entry')->[0]->attributes->getNamedItem('dataset')->nodeValue;
	$info{name} = $entry->getElementsByTagName('name')->[0]->firstChild->nodeValue;
	# NB: the field names don't reflect the real meaning
	if (my $p = $entry->getElementsByTagName('protein')) {
		if ($p->[0]->childNodes and $p->[0]->childNodes->[1]->childNodes) {
			$info{proteinName} = $p->[0]->childNodes->[1]->childNodes->[1]->firstChild->nodeValue;
		}
	}
	if (my $pe_tag = $entry->getElementsByTagName('proteinExistence')) {
		$info{proteinType} = $pe_tag->[0]->attributes->getNamedItem('type')->nodeValue;
	}

	if (my @genes = $entry->getElementsByTagName('gene')) {
		my $name_tags = $genes[0]->getElementsByTagName('name');
		$info{geneName} = $name_tags->[0]->firstChild->nodeValue;
		$info{geneType} = $name_tags->[0]->attributes->getNamedItem('type')->nodeValue;
	} elsif (my @gene_location_list = $entry->getElementsByTagName('geneLocation')) {
		if (my $name_tags = $gene_location_list[0]->getElementsByTagName('name')) {
			if ($name_tags->[0]->firstChild) {
				$info{geneName} = $name_tags->[0]->firstChild->nodeValue;
				if ($info{geneType} = $name_tags->[0]->attributes->getNamedItem('type')) {
					$info{geneType} = $name_tags->[0]->attributes->getNamedItem('type')->nodeValue;
				}
			}
		}
	}

	my @db_refs;
	foreach my $db_ref ($entry->getElementsByTagName('dbReference')) {
		my $db_ref_type = $db_ref->attributes->getNamedItem('type')->nodeValue;
		my $db_ref_id = $db_ref->attributes->getNamedItem('id')->nodeValue;

		$info{dbRefId} = $db_ref_id if $db_ref_type eq 'GeneId';

		my @property_list = $db_ref->getElementsByTagName('property');
		next unless @property_list;
		my $property_type = $property_list[0]->attributes->getNamedItem('type')->nodeValue;
		next if $property_type ne 'entry name';
		
		my $db_ref_name = $property_list[0]->attributes->getNamedItem('value')->nodeValue;

		push(@db_refs, {accession => $info{accession},
						dbRefId => $db_ref_id,
						dbRefType => $db_ref_type,
						dbRefName => $db_ref_name});
	}

	s/\|/\\|/g for values(%info);
	my @l; push(@l, $info{$_}) for qw(accession name dataset proteinName proteinType geneType geneName dbRefId);
	$uniprot_h{$info{accession}} = join('|', @l);

	foreach my $ref (@db_refs) {
		s/\|/\\|/g for values(%$ref);
		my @l; push(@l, $$ref{$_}) for qw(accession dbRefId dbRefType dbRefName);
		$uniprot_evidence_h{$info{accession}} = join('|', @l);
	}
#print "Entry:\n"; print "\t$_\t\"$info{$_}\"\n" for keys %info; print "\tRefs:\n\t\t"; foreach my $r (@db_refs) { print "$_=$$r{$_}; " for keys %$r; print "\n\t\t"; } print "\n";
}
