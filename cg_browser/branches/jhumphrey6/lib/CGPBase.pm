#!/usr/bin/env perl
#package definition#
package CGPBase;
require 5.005;
use Exporter;
our @ISA = "Exporter";
our @methods = qw(quick gene protein sig_peptide tmhmm getorganisms getcontigs details_page);
our %EXPORT_TAGS = (all => [@methods]);
Exporter::export_ok_tags('all');
#/package definition#
use DBI;
use Bio::DB::GFF;
use Bio::DB::GFF::Feature;
#use Template;
#my $tpl = Template->new({INCLUDE_PATH => '/var/www/nbase2010/html/view',}) or  die "template error\n";
use Data::Dumper;
my $output = "";
my $genome='pipeline';
my $db = Bio::DB::GFF->new( -adaptor=>'dbi::mysqlopt',-dsn=>"dbi:mysql:$genome",-user=>'gbrowse2',-pass=>'Brows3r2');

sub query_get_fids($){ # get a list of feature ids
	my $query = shift;
	my $dbi = DBI->connect("DBI:mysql:database=$genome","gbrowse2","Brows3r2");
	my $sth = $dbi->prepare($query);
	$sth->execute();
	my @fids;
	while ( my $fr = $sth->fetchrow_array()){ 
		push(@fids,$fr);
	}
	return \@fids;
}
sub db_raw_sql($){
	my $query = shift;
	my $dbi = DBI->connect("DBI:mysql:database=$genome","gbrowse2","Brows3r2");
	my $sth = $dbi->prepare($query);
	$sth->execute();
}
sub results_factory($$$$){ # takes a feature object and returns a search results object
	my $f = shift; # seqfeature object
	my $term = shift; # fields to display in Match column, OR search term to match
	my $header = shift; # the header to display in the Match column
	my $type = shift;
	my $description="unknown";
	if($f->has_tag('locus_tag')){$description = ($f->get_tag_values('locus_tag'))[0];}
	if($f->has_tag($term)){$description = ($f->get_tag_values($term))[0];}
	else{
		my @tags=$f->get_all_tags();
		foreach my $tag (@tags){
			
			my $value=($f->get_tag_values($tag))[0];
			if($value && $term && $value =~ /$term/i){
				$description=$value;
				last;
			}
		}
	} 
	$description =~ s/^(.{30}).*$/$1.../ if length $description > 30;
	my $row = {
		type => $f->type,
		contig => $f->seq_id,
		source => $f->source(),
		genome => $genome,
		organism => ($f->get_tag_values('organism'))[0],
		locus => $f->display_name,
		start => $f->start,
		end => $f->end,
		length => $f->length,
		feature_id => $f->id,
		match => $description,
		matchtype => $header 
	};
	if(grep(/^$type$/,qw/quick protein gene/)){
		#$$row{'score'}= abs(length($term)-length($description));
		$$row{'score'}=0;#not working :(
	}
	elsif(grep(/^$type$/,qw/sig_peptide/)){
		$$row{'score'}=0;
		$description = sprintf("max C=%s, mean S=%s",($f->get_tag_values('meanS'))[0],($f->get_tag_values('meanS'))[0]);
		my @overlaps = $f->overlapping_features();
		foreach my $subftr (@overlaps){
			if($subftr->has_tag('gene')){
				$description .= sprintf(", %s",($subftr->get_tag_values('gene'))[0]);
				last;
			}
			elsif($subftr->has_tag('product')){
				$description .= sprintf(", %s",($subftr->get_tag_values('product'))[0]);
				if($subftr->type eq 'CDS'){last;}# a CDS product is good enough
			}
		}
		$$row{'match'}=$description;
	}
	elsif(grep(/^$type$/,qw/tmhmm/)){
		$$row{'score'}=0;
		$$row{'length'} = int ($$row{'length'}/3 + 0.5);
		my @overlaps = $f->overlapping_features();
		my @descriptions;
		foreach my $subftr (@overlaps){
			if($subftr->has_tag('gene')){
				push(@descriptions,($subftr->get_tag_values('gene'))[0]);
				last;
			}
			elsif($subftr->has_tag('product')){
				push(@descriptions,($subftr->get_tag_values('product'))[0]);
				if($subftr->type eq 'CDS'){last;}# a CDS product is good enough
			}
		}
		$$row{'match'}=join(', ',@descriptions);
	}
	return $row;
}

sub quick($){
	my $vars = shift;
	my $term = $$vars{'quick'};
	my (@filters,@organisms,@types,@sources);
	my $selected;
	if(defined $$vars{'f_organisms'}){
		@organisms = split('\|',$$vars{'f_organisms'});
		if(0<@organisms){
			foreach (@organisms){
				s/^(.*)$/'$1'/;
			}
			$selected=join(',',@organisms);
			push(@filters,"fid in (select fid from fattribute_to_feature where fattribute_id=(select fattribute_id from fattribute where fattribute_name='organism') and fattribute_value in ($selected))");
		}
	}
	if(defined $$vars{'f_types'}){
		@types = split('\|',$$vars{'f_types'});
		if(!@types){	
			@types=qw/gene CDS misc_feature misc_structure sig_peptide/;
		}
		if(0<@types){
			foreach (@types){
				s/^(.*)$/'$1'/;
			}
			$selected=join(',',@types);
			push(@filters,"fid in (select fid from fdata where ftypeid in (select distinct ftypeid from ftype where fmethod in ($selected)))");
		}
	}
	if(defined $$vars{'f_sources'}){
		@sources = split('\|',$$vars{'f_sources'});
		if(0<@sources){
			foreach (@sources){
				s/^(.*)$/'$1'/;
			}
			$selected=join(',',@sources);
			push(@filters,"fid in (select fid from fdata where ftypeid in (select distinct ftypeid from ftype where fsource in ($selected)))");
		}
	}
	my @results;
	if($term =~ /^\s*$/ && @filters){
		my @type2source;
		foreach my $type (@types){
			$type=~s/'//g;
			if($type eq 'misc_feature'){
				foreach my $source (@sources){
					$source=~s/'//g;
					push(@type2source,"$type:$source");
				}
			}
			else{
				push(@type2source,"$type");
			}
		}
		my @features;
		if(@organisms){
			my @filter_by_organism;
			foreach my $organism (@organisms){
				$organism=~s/'//g;
				push(@features, $db->features(-types=>@type2source,-attributes=>{organism=>$organism}));
			}
		}
		else{
			@features = $db->features(-types=>@type2source);
		}
		foreach my $ftr (@features){
			push(@results, results_factory($ftr,'','Match','tmhmm'));
		}
	}
	else{
		my @termlist = split(/\s+/, $term);
		foreach $term (@termlist){
			my $query = "select fid from fattribute_to_feature where fattribute_value like '\%$term\%' and fattribute_id not in (select fattribute_id from fattribute where fattribute_name like '\%translation\%' or '\%organism\%')";
			if(@filters){$query .= " and " . join(" and ",@filters);}
			my $fids = query_get_fids($query);
			foreach my $fid (@$fids){
				my $f = $db->get_feature_by_id($fid);
				push(@results, results_factory($f,$term,'Match','quick'));
			}
		}
	}
	return \@results;
}
sub gene($){
	my ($vars) = shift;
	my $min=$$vars{'gene_min'};
	my $max=$$vars{'gene_max'};
	my $term=$$vars{'gene_name'};
	$min ||= 0;
	$max ||= 0;
	my @results;
##APPLY SEARCH FILTERS
	my @filters;
	if(defined $$vars{'f_organisms'}){
		my @organisms = split('\|',$$vars{'f_organisms'});
		if(0<@organisms){
			foreach (@organisms){
				s/^(.*)$/'$1'/;
			}
			my $selected=join(',',@organisms);
			push(@filters,"fid in (select fid from feature2tag where type='gene' and tag='organism' and value in ($selected)");
		}
	}
##
	my @termlist = split(/\s+/, $term);
	foreach $term (@termlist){
		my $query = "select fid from feature2tag where type='gene' and tag='gene' and value like '\%$term\%'";
		if(@filters){$query.=" and " . join(" and ",@filters) . ")";}
		my $fids = query_get_fids($query);
		foreach my $fid (@$fids){
			my $f = $db->get_feature_by_id($fid);
			if($min > 0 && $f->length lt $min){next;}
			if($max > 0 && $f->length gt $max){next;}
			push(@results, results_factory($f,$term,"Match",'gene'));
		}
	}
	#my @sorted = sort {$b->{overlap} <=> $a->{overlap}} @results;
	my @sorted = sort {$a->{score} <=> $b->{score}} @results;
	return \@sorted;
}
sub protein($){
	my ($vars) = shift;
	my $min=$$vars{'prot_min'};
	my $max=$$vars{'prot_max'};
	my $term=$$vars{'prot_name'};
	$min ||= 0;
	$max ||= 0;
	my @results;
##APPLY SEARCH FILTERS
	my @filters;
	if(defined $$vars{'f_organisms'}){
		my @organisms = split('\|',$$vars{'f_organisms'});
		if(0<@organisms){
			foreach (@organisms){
				s/^(.*)$/'$1'/;
			}
			my $selected=join(',',@organisms);
			push(@filters,"fid in (select fid from feature2tag where type='CDS' and tag='organism' and value in ($selected)");
		}
	}
##
	my @termlist = split(/\s+/, $term);
	foreach $term (@termlist){
		my $query = "select fid from feature2tag where type='CDS' and tag='product' and value like '\%$term\%'";
		if(@filters){$query.=" and " . join(" and ",@filters) . ")";}
		my $fids = query_get_fids($query);
		foreach my $fid (@$fids){
			my $f = $db->get_feature_by_id($fid);
			if($min > 0 && $f->length lt $min){next;}
			if($max > 0 && $f->length gt $max){next;}
			push(@results, results_factory($f,$term,"Match",'protein'));
		}
	}
	#my @sorted = sort {$b->{overlap} <=> $a->{overlap}} @results;
	my @sorted = sort {$a->{score} <=> $b->{score}} @results;
	return \@sorted;
}
sub sig_peptide($){
	my ($vars) = shift;
	my $meanS=$$vars{'meanS'};
	my $maxC=$$vars{'maxC'};
	my @results;
	my $maxCquery = "select fid from fattribute_to_feature where fattribute_value <= $maxC and fattribute_id=(select fattribute_id from fattribute where fattribute_name='maxC')";
	my $meanSquery = "select fid from fattribute_to_feature where fattribute_value <= $meanS and fattribute_id=(select fattribute_id from fattribute where fattribute_name='meanS')";
	my $query;
	$maxC ||= 0;
	$meanS ||= 0;
	if($maxC && $meanS){$query = $maxCquery . " and fid in ($meanSquery)";}
	elsif($maxC){$query = $maxCquery;}
	elsif($meanS){$query = $meanSquery;}
##APPLY SEARCH FILTERS
	my @filters;
	if(defined $$vars{'f_organisms'}){
		my @organisms = split('\|',$$vars{'f_organisms'});
		if(0<@organisms){
			foreach (@organisms){
				s/^(.*)$/'$1'/;
			}
			my $selected=join(',',@organisms);
			push(@filters,"fid in (select fid from feature2tag where type='sig_peptide' and tag='organism' and value in ($selected)");
		}
	}
	if(@filters){$query.=" and " . join(" and ",@filters) . ")";}
##
	my $fids = query_get_fids($query);
	foreach my $fid (@$fids){ 
		my $f = $db->get_feature_by_id($fid);
		my $searchresult = results_factory($f,'maxC','Match','sig_peptide');
		push(@results, $searchresult);
		
	}
	my @sorted = sort {$b->{length} <=> $a->{length}} @results;
	return \@sorted;
}
sub tmhmm{
	my ($vars) = shift;
	my $minlength=$$vars{'minlength'};
	my $maxlength=$$vars{'maxlength'};
	my $mincount=$$vars{'mincount'};
	my $maxcount=$$vars{'maxcount'};

	my @features=();
	my $lengthquery = '';
	my $countquery = '';
	if($minlength || $maxlength){
		$lengthquery="select fid from feature2tag where type='misc_structure' and tag='length'";
		if($minlength && $minlength>=0){$lengthquery.=" and value>=$minlength"}
		if($maxlength && $maxlength>=0){$lengthquery.=" and value<=$maxlength"}
	}
	if($mincount || $maxcount){
		$countquery="select fid from feature2tag where type='misc_structure' and tag='predicted_number'";
		if($mincount && $mincount>=0){$countquery.=" and value>=$mincount"}
		if($maxcount && $maxcount>=0){$countquery.=" and value<=$maxcount"}
	}
	if($lengthquery){
		$query=$lengthquery;
		if($countquery){$query.=" and fid in ($countquery)";}
	}
	elsif($countquery){
		$query=$countquery;
	}
	else{
		$query="select fid from feature2tag where type='misc_structure'";
	}
	my @results;
##APPLY SEARCH FILTERS
	my @filters;
	if(defined $$vars{'f_organisms'}){
		my @organisms = split('\|',$$vars{'f_organisms'});
		if(0<@organisms){
			foreach (@organisms){
				s/^(.*)$/'$1'/;
			}
			my $selected=join(',',@organisms);
			push(@filters,"fid in (select fid from feature2tag where type='misc_structure' and tag='organism' and value in ($selected)");
		}
	}
	if(@filters){$query.=" and " . join(" and ",@filters) . ")";}
##
	my $fids = query_get_fids($query);
	foreach my $fid (@$fids){
		my $ftr = $db->get_feature_by_id($fid);
		push(@results, results_factory($ftr,'location','feature','tmhmm'));
	}
	return \@results;
}
sub fasta($$){
	my ($source,$fid) = @_;
	my $db = gff_connect($source);
	my $ftr = $db->get_feature_by_id($fid);
	my $fasta=Bio::SeqIO->new(-format=>'fasta');
	my $seq2fasta=$ftr->seq;
	if($ftr->has_tag('locus_tag')){
		my $locus_tag = ($ftr->get_tag_values('locus_tag'))[0];
		my $defline = sprintf("lcl|%s|%s|%d|%d",$locus_tag,$ftr->primary_tag,$ftr->start,$ftr->end);
		$seq2fasta->display_name($defline);
		if($ftr->has_tag('gene')){$seq2fasta->desc(($ftr->get_tag_values('gene'))[0]);}
		elsif($ftr->has_tag('product')){$seq2fasta->desc(($ftr->get_tag_values('product'))[0]);}
		elsif($ftr->primary_tag =~ /gene/){}#skip it! avoids printing duplicate gene/CDS pairs.
		else{$seq2fasta->desc("predicted cds");}
	}
	$fasta->write_seq($seq2fasta) or die "$!\n";
	$seq2fasta->seq($seq2fasta->translate->seq);
	$fasta->write_seq($seq2fasta) or die "$!\n";
	return;
}

sub genome2fasta{
	my $organism=shift;
	my $outfile="$organism.fna";
	my $db = gff_connect('pipeline');
	my $fasta=Bio::SeqIO->new(-file=>">$outfile",-format=>'fasta')or die $!;
	my @features=$db->features(-type=>'CDS',-attributes=>{organism=>$organism});
	foreach my $ftr(@features){
		my $seq2fasta=$ftr->seq;
		my $locus_tag = ($ftr->get_tag_values('locus_tag'))[0];
		my $defline = sprintf("lcl|%s|%s|%d|%d",$locus_tag,$ftr->primary_tag,$ftr->start,$ftr->end);
		$seq2fasta->display_name($defline);
		if($ftr->has_tag('gene')){$seq2fasta->desc(($ftr->get_tag_values('gene'))[0]);}
		elsif($ftr->has_tag('product')){$seq2fasta->desc(($ftr->get_tag_values('product'))[0]);}
		elsif($ftr->primary_tag =~ /gene/){}#skip it! avoids printing duplicate gene/CDS pairs.
		else{$seq2fasta->desc("predicted cds");}
		$fasta->write_seq($seq2fasta) or die "$!\n";
	}
	return $outfile;
}

	
	
sub details_page($$){
	my ($source,$fid) = @_;
	my $db = gff_connect($source);
	my $f = $db->get_feature_by_id($fid);
	my $contig=$f->seq_id;
	my $trans;
	my $seq=$f->seq->seq;
	if($f->strand < 0 ){
		$trans=$f->seq->translate->seq;
	}
	else{
		$trans=$f->seq->translate->seq;
	}
	my @lines = ();
	$trans =~ s/(.{9})(.)/$1$2 /g;
	while ($trans){push(@lines,substr($trans,0,65));if( 66 > length($trans) ){last;}else {$trans=substr($trans,66,-1);}}
	$trans = join("\n",@lines);
	@lines = ();
	$seq =~ s/(.{9})(.)/$1$2 /g;
	while ($seq){push(@lines,substr($seq,0,65));if( 66 > length($seq) ){last;}else {$seq=substr($seq,66,-1);}}
	$seq = join("\n",@lines);
	my $ftr = {start=>$f->start,end=>$f->end,contig=>$contig,source=>$source,sequence=>$seq,translation=>$trans,strand=>$f->strand};
	if ($f->strand <= 0){
		$$ftr{start}=$f->end;
		$$ftr{end}=$f->start;
	}
	my %titles = (gene=>'gene',CDS=>'cds',misc_feature=>'family/domain',sig_peptide=>'signal peptide',misc_structure=>'structure');

	foreach my $ftrtype (keys %titles){
		my @features = $db->overlapping_features(-refseq => $contig, -start => $f->start, -stop => $f->end, -type=>$ftrtype );
		my $i = 0;
		foreach my $subftr (@features){
			my $ftrtags=tags2hash($subftr);
			$$ftrtags{'fid'}=$subftr->id;
			if ($subftr->primary_id() == $f->primary_id()){
				$$ftrtags{'main'}=1;
				if($ftrtype eq 'gene' && defined $$ftrtags{'gene'}){$titles{'gene'}=$$ftrtags{'gene'};}
				if(defined $$ftrtags{'product'}){$titles{$ftrtype}=$$ftrtags{'product'};}
				$$ftr{'title'}=sprintf("%s: %s - %s",$$ftrtags{'organism'},$titles{$ftrtype},$$ftrtags{'locus_tag'});
			}
			if($f->strand <= 0){
				$$ftrtags{'pos'}=$subftr->end . ".." . $subftr->start;
			}
			else{
				$$ftrtags{'pos'}=$subftr->start . ".." . $subftr->end;
			}
			if(!defined($$ftrtags{'length'})){$$ftrtags{'length'}=abs($subftr->end-$subftr->start)+1;}
			if($ftrtype eq 'misc_structure'){
				my @domains=$subftr->get_tag_values('region');
				foreach my $domain (@domains){
					my ($domlocation,$domstart,$domend) = split('\|',$domain);
					push(@{$$ftr{$ftrtype}},{location=>$domlocation,pos=>"$domstart..$domend"});
					$i++;
				}
				next;
			}
			else{
				$$ftr{$ftrtype}[$i]=$ftrtags; 
				$i++;
			}
		}
	}
	return $ftr;
}
sub tags2hash($){
    my $ftr = shift;
    my $data = {};
    my @tags=$ftr->get_all_tags();
    foreach $tag (@tags){ my $tagname = $tag; $tagname =~ s/^\s+//;$tagname =~s/-/_/g;$$data{$tagname}=($ftr->get_tag_values($tag))[0];}
    return $data;
}
sub getorganisms(){
	my $organisms = query_get_fids("select fref from fdata where ftypeid=(select ftypeid from ftype where fmethod='organism')"); 
	return $organisms;
}
sub gettypes(){
	return query_get_fids("select distinct fmethod from ftype");
}
sub getsources(){
	return query_get_fids("select distinct fsource from ftype");
}
sub getcontigs ($){
	my $genome = shift;	
	my @values = ();
	my $db = gff_connect('pipeline');
	my @features = $db->features(-types=>('contig'),-attributes=>{organism=>"$genome"});
	foreach my $feature (@features){
		push(@values,sprintf("%s",$feature->name()));
	}
	@values = sort @values;
	return \@values;
}
sub gff_connect($){
	my $dbname = shift;
	return Bio::DB::GFF->new( -adaptor=>'dbi::mysqlopt',-dsn=>"dbi:mysql:$dbname",-user=>'gbrowse2',-pass=>'Brows3r2') or die "ERROR in CGPBase\n";
}
	
sub blast($){
	use Bio::Tools::Run::StandAloneBlast;
	use IO::String;
	my $vars=shift;#organisms,query(raw sequence;
	my $query=$$vars{'blast'};
	my $seq;
	if($query=~/::::/){$query=join("\n",split("::::",$query));}
	if($query =~ /^\s*>/){#defline - fasta format
		my $strio=IO::String->new($query);
		my $seqio=Bio::SeqIO->new(-fh=>$strio,-format=>'fasta');
		$seq=$seqio->next_seq();
	}
	else{
		$query=~s/[^A-Za-z]//g;
		$seq=Bio::Seq->new(-seq=>$query);
	}
	my @hits;
	my @organisms=split('\|',$$vars{'f_organisms'});
	eval{
		foreach my $organism (@organisms){
	#print STDERR "running blast on $organism evalue $$vars{evalue} program $$vars{program}\n";
			my $db=join("/",$$vars{'blastroot'},$organism);
			my $blast=Bio::Tools::Run::StandAloneBlast->new(database=>$db,program=>$$vars{'program'},e=>$$vars{'evalue'});
			$queries=$blast->blastall($seq);
			while(my $query=$queries->next_result){
				while(my $hit=$query->next_hit){
					while(my $seg=$hit->next_hsp){
						my (undef,$locus_tag)=split('\|',$hit->name);
						push(@hits,{locus_tag=>$locus_tag,desc=>$hit->description,identity=>int(100 * $seg->frac_identical),evalue=>$seg->evalue});
					}
				}
			}
		}
	};
	if( $@ ) {
		$$vars{'message'}=join("\n",$@);
	}
	if(!@hits){
		$$vars{'message'}="Your sequence returned no hits.";
	}
	my @final;
	foreach my $hit(@hits){
		my $locus_tag=$hit->{locus_tag};
		my $query="select fid from feature2tag where tag='locus_tag' and value='$locus_tag' and type='gene'";
		my $fid=query_get_fids($query);
		$hit->{fid}=shift @$fid;
		push (@results,$hit);
	}
	return \@results;
}
sub organism_data{
	my $db = gff_connect('pipeline');
	my @orgs=$db->features(-type=>'organism');
	my %tag_data;
	foreach my $org(@orgs){
		my $tags=tags2hash($org);
		$tag_data{$org->seq_id}=$tags;
	}
	#print STDERR Dumper(%tag_data);
	return \%tag_data;
}		

		


















