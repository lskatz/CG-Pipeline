#!/usr/bin/perl
# usage: [scriptname] prediction.gb *.blast.sql *.ipr_matches.sql vfdb_hits.sql *.signalp_hmm.sql *.signalp_nn.sql *.tmhmm_location.sql *.tmhmm.sql > out.gbk
# also writes annotation.gff

use warnings;
use Bio::Seq;
use Bio::Seq::RichSeq;
use Bio::SeqIO;
use Bio::Tools::GFF;
use File::Basename;
use Data::Dumper;
sub gene_factory($$$$){
	my ($start,$end,$strand,$tags) = @_;
#		return Bio::SeqFeature::Generic->new(-primary=>'gene',-start=>$start,-end=>$end,-strand=>$strand,-tag=>{locus_tag=>$locus_tag,gene=>$gene});
	return Bio::SeqFeature::Generic->new(-primary=>'gene',-start=>$start,-end=>$end,-strand=>$strand,-tag=>$tags);
}
sub remove_tags($){
	my $ftr = shift;
	my @note=();
	my @tags=$ftr->get_all_tags();
	foreach $tag (@tags){
		@note = (@note,
				("$tag:" . ($ftr->get_tag_values($tag))[0])
			);
	}
	foreach $tag (@tags){
		$ftr->remove_tag($tag);
	}
	return join("; ",@note);
}
sub replace_tags($){
	my $ftr = shift;
	if(!$ftr->has_tag('note')){return;}
	my @note = split(';',join('; ',$ftr->get_tag_values('note')));
	foreach (@note){
		my ($tag,$value)=split(':');
		$ftr->add_tag_value($tag,$value);
	}
	$ftr->remove_tag('note');
}
sub start2nuc($$){
	my ($start,$coord) = @_;
	return $start - 3 + 3*$coord;
}
sub end2nuc($$){
	my ($start,$coord) = @_;
	return $start - 1 + 3*$coord;
}
sub tags2hash($){
	my $ftr = shift;
	my %data = ();
	my @tags=$ftr->get_all_tags();
	foreach $tag (@tags){$data{$tag}=($ftr->get_tag_values($tag))[0];}
	return \%data;
}
sub addtranslation($$){
	my $ftr = shift;
	my $seq = shift;
	if($ftr->strand < 0){
		$ftr->add_tag_value('translation',Bio::Seq->new(-seq=>$seq->subseq($ftr->start,$ftr->end))->revcom->translate->seq);
	}
	else{
		$ftr->add_tag_value('translation',Bio::Seq->new(-seq=>$seq->subseq($ftr->start,$ftr->end))->translate->seq);
	}
}
	
	
my ($prediction,@sqlfiles)=@ARGV;
my $gbin = Bio::SeqIO->new(-file=>"$prediction",-format=>'genbank');
my $gbout = Bio::SeqIO->new(-format=>'genbank');
my $gff = Bio::Tools::GFF->new(-file=>">annotation.gff",-gff_version => 3);
my %data;
my $locus_tag;
my %spfeats;
my @args;
my $newftr;
my $type;
my @params = ();
my $tag;

my %map = (
	blast => [qw/locus_tag uniprot_id length name rank score bits evalue identity positives/],
	ipr_matches => [qw/locus_tag accession_num product database_name start end evalue status evidence/],
	signalp_hmm	=> [ qw/locus_tag prediction signal_peptide_probability max_cleavage_site_probability start end/],
	signalp_nn => [ qw/locus_tag measure_type start end value cutoff is_signal_peptide/],
	tmhmm_location => [qw/locus_tag location start end/],
	tmhmm => [ qw/locus_tag predicted_number length expected_number_aa expected_number_aa_60 total_prob_n_in/],
	#ipr_hits => [qw/locus_tag length domain_id product type/],
	#ipr_classifications => [qw/locus_tag go_id class_type category description/],
	#ipr_childrefs => [qw/locus_tag interpro_hit_id interpro_child_reference/],
	vfdb_hits => [qw/locus_tag target_id evalue coverage db_name identity/]
);
STAGE_1:
#Stage 1: Make a hash of loci to SeqFeature objects from each data line
foreach $sqlfile(@sqlfiles){
	$type = (reverse split('\.', basename($sqlfile,'.sql')))[0];
	if(!defined($type) || !defined($map{$type})){ print STDERR "Input filename rejected: $type\n";next;}
	@params = @{$map{$type}};
	open (FH, "<$sqlfile") or die $!; #one prediction file containing locus tags
	while (<FH>){
		chomp;
		# skip a row of empties
		if( /^\s*\|/ ){ print STDERR "Skipping blank line:$_\n"; next;}
		$newftr = Bio::SeqFeature::Generic->new(-primary=>$type,-start=>'0',-end=>'0');
		s/\\\|/::/g;
		@values = split('\|');
		if(scalar @values != scalar @params){print STDERR "Wrong number of fields:$_\n";next;}
		for(my $i=0;$i<scalar @params;$i++){
			if($params[$i] eq "start"){
				$newftr->start($values[$i]);
			}
			elsif($params[$i] eq "end"){
				$newftr->end($values[$i]);
			}
			else{
				#if($params[$i] eq "locus_tag"){
				#	$values[$i] =~ s/_Pipeline_draft_/_/g;
				#}
				$newftr->add_tag_value($params[$i],$values[$i]);
			}
		}
		if(!defined($spfeats{$values[0]})){$spfeats{$values[0]} = [];}
		push (@{$spfeats{$values[0]}}, $newftr );
	}
	close(FH);
}
STAGE_2:
#Stage 2: Select a name for each locus from available SeqFeature objects and integrate the other subfeatures
while( my $seq=$gbin->next_seq()){ #each contig
	my $strain=$seq->primary_seq->desc;
	$strain =~ s/^([^_]+).*$/$1/;
	my $contig=$seq->primary_seq->display_id;
	$contig =~ s/contig[0]*([0-9]+)_.*$/$1/;
	$seq->display_id(sprintf("%s_%04d",$strain,$contig));
	my @feats = $seq->all_SeqFeatures(); # features here are loci in the genome
	$seq->flush_SeqFeatures();
	foreach my $ftr (@feats){ # each feature in the contig
		#print $ftr->primary_tag() . "\n";next;
		if($ftr->primary_tag() eq 'gene'){next;} # we get two features, gene and CDS, duplicates
		$locus_tag = ($ftr->each_tag_value('locus_tag'))[0];
		#$locus_tag =~ s/_Pipeline_draft_/_/g;
		#print STDERR $ftr->gff_string($gff) . "\n";
		if (!(defined($spfeats{$locus_tag}) && (0<scalar @{$spfeats{$locus_tag}}))){next;}
		#initialize the naming hash
		%data = (type=>'none',name=>'hypothetical',product=>'putative',identity=>1,evalue=>999,count=>'0');
		foreach $newftr ( @{$spfeats{$locus_tag}}){
			#if($data{'product'} !~ /hypothetical|putative|probable|conserved/i){
			#	#the gene has been satisfactorily named by a uniprot blast hit
			#	next;	
			#	#otherwise, the gene name will be the concatenation of all other evidence
			#}
			$newftr->strand($ftr->strand);
#User should always list data files in order: blast,ipr_matches,signalp,tmhmm
			if($newftr->primary_tag eq 'blast'){
				$newftr->start($ftr->start());
				$newftr->end($ftr->end());
				$newftr->strand($ftr->strand());
				my $identity=($newftr->get_tag_values('identity'))[0];
				my $evalue=($newftr->get_tag_values('evalue'))[0];
				my $product=($newftr->get_tag_values('name'))[0];
				my $source =$newftr->source_tag();
				my $locus_tag=($newftr->get_tag_values('locus_tag'))[0];
				my $gene="unnamed";
				if($identity<91 || $evalue>1e-9){next;}
				if(	($data{'evalue'} >= $evalue || #compare this blast hit to the previous one
					$data{'identity'} <= $identity ) ){
					%data=(type=>'blast',product=>$product,identity=>$identity,evalue=>$evalue,locus_tag=>$locus_tag,count=>$data{'count'}++);
					if($product =~ /GN=([^=]+)\s+[A-Z]{2,}=/){$gene=$1;}# look for a short name
					if($product =~ /^([^=]+)\s+[\S^=]+=.*$/){$data{'product'}=$1;}
					my $note=remove_tags($newftr);
					$newftr->add_tag_value('note',$note);
					$newftr->add_tag_value('product',$data{'product'});
					$newftr->add_tag_value('locus_tag',$data{'locus_tag'});
					$newftr->primary_tag('CDS');
					addtranslation($newftr,$seq);
					$newftr->display_name($data{'name'});
#remove any existing feature for this locus
					my @currentfeatures=$seq->remove_SeqFeatures();
					foreach my $curftr(@currentfeatures){
						if(!$curftr->has_tag('locus_tag')){next;}
						my $protid=($curftr->get_tag_values('locus_tag'))[0];
						if($locus_tag !~ /^$protid$/){
							$seq->add_SeqFeature($curftr);
						}
					}
#create the main feature
					my $parent = gene_factory($newftr->start,$newftr->end,$newftr->strand,{locus_tag=>$locus_tag});
					$parent->add_SeqFeature($newftr);
					if($gene ne 'unnamed' && $gene !~ /[_-]+/ && $gene !~ /^[A-Z]+/ && $gene !~ /^[0-9]{2,}/){$parent->add_tag_value('gene',$gene);}
					$seq->add_SeqFeature($parent);
					$ftr=$parent;#future features (below, from ipr_matches, etc) will be subfeatures of this one
				}
			}
			elsif($newftr->primary_tag eq 'ipr_matches'){
				my $product=($newftr->get_tag_values('product'))[0];
				my $evidence=($newftr->get_tag_values('evidence'))[0];
				my $note = remove_tags($newftr);
				$newftr->add_tag_value('note',$note);
				$newftr->add_tag_value('product',$product);
				$newftr->primary_tag('misc_feature');
				$newftr->start(start2nuc($ftr->start,$newftr->start));
				$newftr->end(end2nuc($ftr->start,$newftr->end));
				$ftr->add_SeqFeature($newftr);
			}
			elsif($newftr->primary_tag eq 'vfdb_hits'){
				my @currentfeatures = $ftr->get_SeqFeatures();
				foreach my $each_ftr (@currentfeatures){
					if ($each_ftr->primary_tag() eq 'CDS'){
						$each_ftr->add_tag_value('vfdb_id',($newftr->get_tag_values('target_id'))[0]);
					}
				}
			}
			elsif($newftr->primary_tag eq 'signalp_nn'){
			#	if(($newftr->remove_tag('is_signal_peptide'))[0] eq 'NO'){next;}
				my $parent;
				my $aa_pos = $newftr->end ?  sprintf("%d-%d",$newftr->start, $newftr->end): sprintf("%d",$newftr->start);
				$newftr->start(start2nuc($ftr->start,$newftr->start));
				if($newftr->end){$newftr->end(end2nuc($ftr->start,$newftr->end));}
				$newftr->remove_tag('locus_tag');
				my $tags = tags2hash($newftr);
				my @currentfeatures = $ftr->get_SeqFeatures();
				foreach my $each_ftr (@currentfeatures){
					if ($each_ftr->primary_tag() eq 'sig_peptide'){#there can be only one
						$parent=$each_ftr;
						last;
					}
				}
				if(!defined($parent)){
					$parent=Bio::SeqFeature::Generic->new(-primary=>'sig_peptide',-start=>$ftr->start,-end=>$ftr->end,-strand=>$ftr->strand,tag=>{locus_tag=>$locus_tag});
					$ftr->add_SeqFeature($parent);
				}
				#add tags to parent as appropriate
				#/maxC=[pos],value=[value],cutoff=[cutoff]
				#/maxY=[pos],value=[value],cutoff=[cutoff]
				#/maxS=[pos],value=[value],cutoff=[cutoff]
				#/meanS=[pos],value=[value],cutoff=[cutoff]
				#D:give the coordinates to parent
				my $measure_tag = $$tags{'measure_type'};
				$measure_tag =~ s/[\. ]//g;
				if($measure_tag eq 'D'){
					if($$tags{'is_signal_peptide'} eq 'NO'){
					#remove this feature from the locus, it has no signal peptide
						@currentfeatures = $ftr->remove_SeqFeatures();
						foreach my $each_ftr (@currentfeatures){
							if ($each_ftr->primary_tag() eq 'sig_peptide'){next;}
							else {$ftr->add_SeqFeature($each_ftr);}
						}
						next;
					}
					$parent->start($newftr->start);
					$parent->end($newftr->end);
					#addtranslation($parent);
					$measure_tag='D-score';
				}
				$parent->add_tag_value($measure_tag,sprintf("aa_pos %s,value %s,cutoff %s",$aa_pos,$$tags{'value'},$$tags{'cutoff'}));
			}
					
			elsif($newftr->primary_tag eq 'signalp_hmm'){
				if($newftr->start < 1){next;}
				my $parent;
				$newftr->remove_tag('locus_tag');
				my $tags = tags2hash($newftr);
				my $aa_pos = sprintf("%d-%d",$newftr->start,$newftr->end);
				$newftr->start(start2nuc($ftr->start,$newftr->start));
				$newftr->end(end2nuc($ftr->start,$newftr->end));
				my @currentfeatures = $ftr->get_SeqFeatures();
				foreach my $each_ftr (@currentfeatures){
					if(!$each_ftr->has_tag('locus_tag')){next;}
					if (($each_ftr->get_tag_values('locus_tag'))[0] eq $locus_tag && $each_ftr->primary_tag() eq 'sig_peptide'){
						$parent=$each_ftr;
						last;
					}
				}
				if(!defined($parent)){
					$parent=Bio::SeqFeature::Generic->new(-primary=>'sig_peptide',-start=>$ftr->start,-end=>$ftr->end,-strand=>$ftr->strand,tag=>{locus_tag=>$locus_tag});
					$ftr->add_SeqFeature($parent);
				}
				#hmm:
				#/cleavage_site=[start-end]
				#/cleavage_prob=[cleav site prob]
				#/sigp_prob=[prob]
				#and give parent hmm's prediction
				$parent->add_tag_value('cleavage_site',$aa_pos);
				$parent->add_tag_value('cleavage_prob',$$tags{'max_cleavage_site_probability'});
				$parent->add_tag_value('signalp_prob',$$tags{'signal_peptide_probability'});
				$parent->add_tag_value('product',$$tags{'prediction'});
			}
			elsif($newftr->primary_tag eq 'tmhmm_location'){
				$newftr->remove_tag('locus_tag');
				my $note = remove_tags($newftr);
				$newftr->add_tag_value('note',$note);
				$newftr->add_tag_value('locus_tag',$locus_tag);
				$newftr->primary_tag('misc_structure');
				$newftr->start(start2nuc($ftr->start,$newftr->start));
				$newftr->end(end2nuc($ftr->start,$newftr->end));
				#addtranslation($newftr);
				$ftr->add_SeqFeature($newftr);
			}
			elsif($newftr->primary_tag eq 'tmhmm'){
				my @currentfeatures = $ftr->get_SeqFeatures();
				foreach my $each_ftr (@currentfeatures){
					if ($each_ftr->primary_tag() eq 'misc_structure'){
						$newftr->remove_tag('locus_tag');
						my $note = ($each_ftr->get_tag_values('note'))[0] . "; " . remove_tags($newftr); 
						$each_ftr->remove_tag('note');
						$each_ftr->add_tag_value('note',$note);
						$newftr->add_tag_value('locus_tag',$locus_tag);
						last;
					}
				}
			}
		}
#post-processing
		my @features = sort {$a->start<=>$b->start}$seq->all_SeqFeatures();
		my @finalfeatures = ();
		foreach $newftr(@features){
			if($newftr->has_tag('gene')){ 
				my $gene = ($newftr->get_tag_values('gene'))[0];
				if($gene eq 'unnamed' ||
					$gene =~ /^[A-Z]+/ ||
					$gene =~ /^[0-9]{2,}/){
					$newftr->remove_tag('gene');
				}
			}
			push(@finalfeatures,$newftr);
		}
		$seq->flush_SeqFeatures();
		$seq->add_SeqFeature(@finalfeatures);
	}
	$gbout->write_seq($seq);
	foreach my $ftr ($seq->all_SeqFeatures()){
		$ftr->seq_id($seq->display_id());
		replace_tags($ftr);
		$gff->write_feature($ftr);
	}
}
