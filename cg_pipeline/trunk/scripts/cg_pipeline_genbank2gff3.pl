#!/usr/bin/env perl
# generate fasta file for an organism containing organism feature, contigs, and predicted loci

use strict;
use warnings;
use Data::Dumper;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::GFF;
use File::Basename;
sub note2tags($);
my %type2source=(
	gene=>"UniProt",
	CDS=>"UniProt",
	misc_feature=>"InterProScan",
	misc_structure=>"TMHMM",
	sig_peptide=>"SignalP"
);
	
if(!@ARGV){print "$0 annotation.gb\n";exit;}
	
my ($annotation) = @ARGV;
my $ftr = Bio::SeqFeature::Generic->new();

my $gbin = Bio::SeqIO->new(-file=>"<$annotation",-format=>'genbank');
#my $assout = Bio::SeqIO->new(-format=>'fasta');
my $gff = Bio::Tools::GFF->new(-gff_version => 3);
my $fasta = Bio::SeqIO->new(-format=>'fasta');

# 1: contigs
my $org;
while( my $gbseq=$gbin->next_seq()){ #each contig
	if(!$org){
		$org=$gbseq->primary_seq->desc;
		if($org =~ /^([^,]+),/){$org =$1;}
		$org =~ s/_Pipeline.*$//;
		$gff->write_feature(Bio::SeqFeature::Generic->new(-seq_id=>$org,-source=>'CGPipeline',-primary_tag=>'organism',-start=>1,-end=>1,tag=>{genus=>$gbseq->species->genus,species=>$gbseq->species->species}));
	}
	my $contig=$gbseq->primary_seq->display_id;
	my $seq = $gbseq->primary_seq;
	$gff->write_feature(Bio::SeqFeature::Generic->new(-seq_id=>$seq->display_id,-source=>'Assembly',-primary_tag=>'contig',-start=>1,-end=>$seq->length,-tag=>{Name=>$seq->display_id,organism=>$org}));
}
$gbin->close;
# 2: features
$gbin = Bio::SeqIO->new(-file=>"<$annotation",-format=>'genbank');
while( my $gbseq=$gbin->next_seq()){ #each contig
	my $contig=$gbseq->primary_seq->display_id;
	my $seq = $gbseq->primary_seq;
	my @features=$gbseq->get_SeqFeatures();
	foreach my $ftr(@features){
		note2tags($ftr);
		my $type=$ftr->primary_tag;
		$ftr->add_tag_value("organism",$org);
		my $source="unknown";
		if($type =~ /^misc_feature$/){
			if($ftr->has_tag("database_name")){
				($source)=$ftr->get_tag_values("database_name");
			}
			elsif($ftr->has_tag("evidence")){
				($source)=$ftr->get_tag_values("evidence");
			}
		}
		else{
			$source=$type2source{$type};
		}
		$ftr->source_tag($source);
		$gff->write_feature($ftr);
	}
}
$gbin->close;
# 3. fasta
print "##fasta\n";
$gbin = Bio::SeqIO->new(-file=>"<$annotation",-format=>'genbank');
while( my $gbseq=$gbin->next_seq()){ #each contig
	my $contig=$gbseq->primary_seq->display_id;
	my $seq = $gbseq->primary_seq;
	my $residues=$seq->seq;
	$residues=~s/N/-/g;
	$seq->seq($residues);
	$seq->display_id($contig);
	$seq->desc('');
	$fasta->write_seq($seq);
}
exit;

sub note2tags($){
	my $ftr = shift;
	if(!$ftr->has_tag('note')){return;}
	my @notes = $ftr->get_tag_values('note');
	$ftr->remove_tag('note');
	foreach my $note (@notes){
		my @attrs=split(/;/,$note);
		foreach my $attr(@attrs){
			my ($tag,$value)=split(/\s*:\s*/,$attr);
			if(! $tag  or ! $value){next;}
			$tag =~ s/^\s+//;
			$tag =~ s/\s+$//;
			$value =~ s/^\s+//;
			$value =~ s/\s+$//;
			$ftr->add_tag_value($tag,$value);
		}
	}
	my @tags = $ftr->all_tags();
	#fix bad tag names
	foreach my $tag (@tags){
		my @values = $ftr->get_tag_values($tag);
		$ftr->remove_tag($tag);
		if(0>=scalar @values){next;}
		$tag =~ s/^\s+//;
		$tag =~ s/\s+$//;
		#remove duplicate values
		my @uniquevalues=("uninitialized");
		foreach my $value(sort @values){
			if($value ne $uniquevalues[$#uniquevalues]){push(@uniquevalues,$value);}
		}	
		shift @uniquevalues;
		if(0<scalar @uniquevalues){
			$ftr->add_tag_value($tag,@uniquevalues);
		}
	}
}
