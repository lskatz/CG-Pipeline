#!/usr/bin/env perl
use warnings;
use strict;
use Bio::Seq;
use Bio::Seq::RichSeq;
use Bio::SeqIO;
use Bio::Tools::GFF;
use File::Basename;
use Data::Dumper;
use Getopt::Long;
	
exit(&main());

sub main(){
	if(@ARGV < 1){
		print STDERR "Usage: ". basename($0) . " --prediction=prediction.gb --inputdir=annotation-sql-folder --gb=genbank-output-file --gff=gff-output-file --fasta=fasta-output-file [--organism=organism_id] \n";
		return 0;
	}
	my %data;
	my %spfeats;
	my @args;
	my $newftr;
	my $type;
	my @params = ();
	my $tag;
	
	my %map = (
		blast => [qw/locus_tag uniprot_id length name rank score bits evalue identity positives/],
		vfdb_hits => [qw/locus_tag target_id evalue coverage db_name identity/],
		ipr_matches => [qw/locus_tag accession_num product database_name start end evalue status evidence/],
		signalp_hmm	=> [ qw/locus_tag prediction signal_peptide_probability max_cleavage_site_probability start end/],
		signalp_nn => [ qw/locus_tag measure_type start end value cutoff is_signal_peptide/],
		tmhmm => [ qw/locus_tag length predicted_number expected_number_aa expected_number_aa_60 total_prob_n_in/],
		tmhmm_location => [qw/locus_tag location start end/],
		uniprot => [qw/ac ac2 source product pro_type gene_type gene_name gene_id/]
		#ipr_hits => [qw/locus_tag length domain_id product type/],
		#ipr_classifications => [qw/locus_tag go_id class_type category description/],
		#ipr_childrefs => [qw/locus_tag interpro_hit_id interpro_child_reference/],
	);

	my $settings = {};
	GetOptions($settings,('organism=s','prediction=s','inputdir=s','gb=s','gff=s','fasta=s')) or die;
	my $organism="organism";
	if(defined($$settings{'organism'})){
		$organism=$$settings{'organism'};
	}
	my @sqlfiles;
	my $atndir;
	opendir($atndir,$$settings{'inputdir'}) or die "unable to open directory $$settings{'inputdir'}:$!\n";
	my @atnfiles=readdir($atndir);
	foreach $type ((qw/uniprot blast ipr_matches signalp_hmm signalp_nn tmhmm.sql tmhmm_location vfdb_hits/)){
		my @files=grep(/$type/,@atnfiles);
		if(@files){push(@sqlfiles,$$settings{'inputdir'} . "/" . $files[0]);}
	}
	
	my $gbin = Bio::SeqIO->new(-file=>$$settings{'prediction'},-format=>'genbank');
	my $gbout = Bio::SeqIO->new(-file=>">".$$settings{'gb'},-format=>'genbank');
	my $fasta;
	if(defined($$settings{'fasta'})){
		$fasta = Bio::SeqIO->new(-file=>">".$$settings{'fasta'},-format=>'');
	}
	my $gff;
	if(defined($$settings{'gff'})){
		$gff = Bio::Tools::GFF->new(-file=>">".$$settings{'gff'},-gff_version => 3);
	}

	my %uniprot;

	STAGE_1:
	#Stage 1: Make a hash of loci to SeqFeature objects from each data line
	my $blastcount=0;
	foreach my $sqlfile(@sqlfiles){
		$type = (reverse split('\.', basename($sqlfile,'.sql')))[0];
		if(!defined($type) || !defined($map{$type})){ print STDERR "Input filename rejected: $type\n";next;}
		@params = @{$map{$type}};
		open (FH, "<$sqlfile") or die $!; #one prediction file containing locus tags
print STDERR "parsing $sqlfile\n";
		while (<FH>){
			chomp;
			# skip a row of empties
			if( /^\s*\|/ ){ print STDERR "Skipping blank line:$_\n"; next;}
			
				
			$newftr = Bio::SeqFeature::Generic->new(-primary=>$type,-start=>'0',-end=>'0');
			s/\\\|/::/g;
			my @values = split('\|');
			if(scalar @values != scalar @params){if($type ne 'uniprot'){print STDERR "Wrong number of fields:" . scalar @values . "\n";next;}}
			if($type eq 'uniprot'){
				@{$uniprot{$values[0]}}{@params}=@values;#Wow! So that's how you make a hash out of two arrays.
				next;
			}
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
		if($organism eq 'organism'){#not set
			$organism=$seq->primary_seq->desc;
			$organism =~ s/^([^_]+).*$/$1/;
		}
		my $contig=$seq->primary_seq->display_id;
		$contig =~ s/[a-zA-Z]*[0]*([0-9]+)_.*$/$1/;
		$seq->display_id(sprintf("%s_%04d",$organism,$contig));
		my @feats = $seq->all_SeqFeatures(); # features here are loci in the genome
		#$seq->flush_SeqFeatures(); # try removing this.....
		foreach my $ftr (@feats){ # each feature in the contig---each locus
			if($ftr->primary_tag() eq 'gene'){next;} # we get two features, gene and CDS, duplicates
			my $locus_tag = ($ftr->get_tag_values('locus_tag'))[0];
			#if (!(defined($spfeats{$locus_tag}) && (0<scalar @{$spfeats{$locus_tag}}))){print STDERR "skipping due to having no features:$locus_tag\n"; next;}
			if (!(defined($spfeats{$locus_tag}))){
#				print STDERR "skipping due to having no features:$locus_tag\n";
				$seq->add_SeqFeature(@feats);# put predictions back in
				next;
			}
			#initialize the naming hash
			%data = (type=>'none',name=>'hypothetical',product=>'putative',identity=>1,evalue=>999,count=>'0');
			my $ftrcount = scalar @{$spfeats{$locus_tag}};
			#print STDERR "$ftrcount features for $locus_tag\n";
			foreach $newftr ( @{$spfeats{$locus_tag}}){
				$newftr->strand($ftr->strand);
				if($newftr->primary_tag eq 'blast'){
					$newftr->start($ftr->start());
					$newftr->end($ftr->end());
					$newftr->strand($ftr->strand());
					my $tags=tags2hash($newftr);
				#	my $identity=($newftr->get_tag_values('identity'))[0];
				#	my $evalue=($newftr->get_tag_values('evalue'))[0];
				#	my $product=($newftr->get_tag_values('name'))[0];
				#	my $locus_tag=($newftr->get_tag_values('locus_tag'))[0];
					$newftr->remove_tag('name');#particularly shitty tag riddled with '='
					my $source =$newftr->source_tag();
					my $gene="unnamed";
					if($$tags{'identity'}<91 || $$tags{'evalue'}>1e-9){
# considering writing the empty features back to the locus and calling this the last iteration
						next;
					}
					if(	($data{'evalue'} >= $$tags{'evalue'} || $data{'identity'} <= $$tags{'identity'} ) ){ #compare this blast hit to the previous one
						$$tags{'product'}=$uniprot{$$tags{'uniprot_id'}}{'product'};
						%data=(type=>'blast',count=>$data{'count'}+1,%$tags);
						$gene=$uniprot{$$tags{'uniprot_id'}}{'gene_name'};
#print STDERR "$gene\n";
						my $note=remove_tags($newftr);
						$newftr->add_tag_value('note',$note);
						$newftr->add_tag_value('product',$data{'product'});
						$newftr->add_tag_value('locus_tag',$data{'locus_tag'});
						$newftr->source_tag('UniProt');
						$newftr->primary_tag('CDS');
						addtranslation($newftr,$seq);
						$newftr->display_name($data{'name'});
	#remove any existing feature for this locus
						my @currentfeatures=$seq->remove_SeqFeatures();
						foreach my $curftr(@currentfeatures){
							if(!$curftr->has_tag('locus_tag')){next;}
							my $protid=($curftr->get_tag_values('locus_tag'))[0];
							if($$tags{'locus_tag'} !~ /^$protid$/){
								$seq->add_SeqFeature($curftr);
							}
						}
	#create the main feature
						my $parent = gene_factory($newftr->start,$newftr->end,$newftr->strand,{locus_tag=>$$tags{'locus_tag'}});
						$parent->add_SeqFeature($newftr);
						if(defined($gene)){
							if(!($gene eq 'unnamed' || ($gene =~ /[_-]+/) || ($gene =~ /^[A-Z]+/) || ($gene =~ /^[0-9]{2,}/))){
								$parent->add_tag_value('gene',$gene);
							}
						}
						$parent->add_tag_value('gene',$gene);
						$seq->add_SeqFeature($parent);
						$ftr=$parent;#future features (below, from ipr_matches, etc) will be subfeatures of this one
					}
				}
				elsif($newftr->primary_tag eq 'ipr_matches'){
					my $tags=tags2hash($newftr);
					my $note = remove_tags($newftr);
					$newftr->add_tag_value('note',$note);
					$newftr->add_tag_value('product',$$tags{'product'});
					$newftr->primary_tag('misc_feature');
					$newftr->source_tag($$tags{'database_name'});
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
						$parent=Bio::SeqFeature::Generic->new(-primary=>'sig_peptide',-source_tag=>'SignalP',-start=>$ftr->start,-end=>$ftr->end,-strand=>$ftr->strand,tag=>{locus_tag=>$locus_tag});
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
						$parent=Bio::SeqFeature::Generic->new(-primary=>'sig_peptide',-source_tag=>'SignalP',-start=>$ftr->start,-end=>$ftr->end,-strand=>$ftr->strand,tag=>{locus_tag=>$locus_tag});
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
				elsif($newftr->primary_tag eq 'tmhmm'){
					$newftr->remove_tag('locus_tag');
					my $note = remove_tags($newftr);
					$newftr->add_tag_value('note',$note);
					$newftr->add_tag_value('locus_tag',$locus_tag);
					$newftr->add_tag_value('product','transmembrane structure');
					$newftr->primary_tag('misc_structure');
					$newftr->source_tag('TMHMM');
					#$newftr->start(start2nuc($ftr->start,$newftr->start));
					#$newftr->end(end2nuc($ftr->start,$newftr->end));
					$newftr->start($ftr->start);
					$newftr->end($ftr->start);# there are no regions yet, added from tmhmm_location.sql
					#addtranslation($newftr);
					$ftr->add_SeqFeature($newftr,'EXPAND');
				}
				elsif($newftr->primary_tag eq 'tmhmm_location'){
					my @currentfeatures = $ftr->get_SeqFeatures();
					foreach my $each_ftr (@currentfeatures){
						if ($each_ftr->primary_tag() eq 'misc_structure'){
							my $tags=tags2hash($newftr);
							$each_ftr->add_tag_value('region',sprintf("%s|%d|%d",$$tags{'location'},$newftr->start,$newftr->end));
							my ($regionstart,$regionend)=(start2nuc($each_ftr->start,$newftr->start),end2nuc($each_ftr->start,$newftr->end));
							if($each_ftr->start == $each_ftr->end){#start expanding the feature
								$each_ftr->start($regionstart);
								$each_ftr->end($regionend);
							}
							else{
								if($each_ftr->start>$regionstart){$each_ftr->start($regionstart);}
								if($each_ftr->end<$regionend){$each_ftr->end($regionend);}
							}
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
					if(defined($gene)){
						if( 
						$gene eq 'unnamed' ||
						$gene =~ /^[A-Z]+/ ||
						$gene =~ /^[0-9]{2,}/){
							$newftr->remove_tag('gene');
						}
					}
				}
				push(@finalfeatures,$newftr);
			}
			$seq->flush_SeqFeatures();
			$seq->add_SeqFeature(@finalfeatures);
		}
		$gbout->write_seq($seq);
# write gff and fasta
		if($gff or $fasta){
			foreach my $ftr ($seq->all_SeqFeatures()){
				$ftr->seq_id($seq->display_name());
				replace_tags($ftr);
				$ftr->add_tag_value('organism',$organism);
				if($gff){
					$gff->write_feature($ftr);
				}
				if($fasta){
					my $seq2fasta = $ftr->seq;
					my $seqid = scalar $seq2fasta->display_name;
					if($ftr->has_tag('locus_tag')){
						my $locus_tag = ($ftr->get_tag_values('locus_tag'))[0];
						my $defline = sprintf("lcl|%s|%s|%d|%d",$locus_tag,$ftr->primary_tag,$ftr->start,$ftr->end);
						$seq2fasta->display_name($defline);
						if($ftr->has_tag('gene')){$seq2fasta->desc(($ftr->get_tag_values('gene'))[0]);}
						elsif($ftr->has_tag('product')){$seq2fasta->desc(($ftr->get_tag_values('product'))[0]);}
						elsif($ftr->primary_tag =~ /gene/){next;}#skip it! avoids printing duplicate gene/CDS pairs.
						else{$seq2fasta->desc("predicted cds");}
						$fasta->write_seq($seq2fasta) or die "$!\n";
					}
				}
			}
		}
	}
	return 0;
}#end main

sub gene_factory($$$$){
	my ($start,$end,$strand,$tags) = @_;
	return Bio::SeqFeature::Generic->new(-primary=>'gene',-source_tag=>'UniProt',-start=>$start,-end=>$end,-strand=>$strand,-tag=>$tags);
}
sub remove_tags($){
	my $ftr = shift;
	my @note=();
	my @tags=$ftr->get_all_tags();
	my $value;
	foreach my $tag (@tags){
		$value=($ftr->get_tag_values($tag))[0];
		$value=~s/=/:/g;
		@note = (@note,
				("$tag:" . $value)
			);
	}
	foreach my $tag (@tags){
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
		if($ftr->has_tag($tag)){next;}#things like 'product'
		else{$ftr->add_tag_value($tag,$value);}
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
	foreach my $tag (@tags){$data{$tag}=($ftr->get_tag_values($tag))[0];}
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
