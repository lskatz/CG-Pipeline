#=============================================================
# $URL$ 
# $Id: CGPBase.pm 36 2010-06-01 14:22:03Z jhumphrey6 $
#=============================================================

package CGPBase;

use 5.006;
use Exporter;
our @ISA = "Exporter";
our @methods = qw(pubmed kegg interpro swissprot keggpathway protein details getcontigs);
our %EXPORT_TAGS = (all => [@methods]);
Exporter::export_ok_tags('all');
#/package definition#
use DBI;
use Bio::DB::GFF;
use Bio::DB::GFF::Feature;
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

sub results_factory($$$){ # takes a feature object and returns a search results object
	my $f = shift; # seqfeature object
	my $query = shift; # the search term you sent to match some tag
	my $header = shift; # the header to display in the Match column
	my @tags = $f->get_all_tags();
	foreach my $tag (@tags){
		my @tagvalues = $f->get_tag_values($tag);
		foreach my $value (@tagvalues)
		{
			if($value =~ /$query/i) # this is the tag matching your search term
			{
				$value =~ s/_/ /g;
				my @notes = $f->get_tag_values($tag);
				my $overlap=scalar( $f->overlapping_features() );
				my $row = {
					overlap => $overlap,
					contig => $f->seq_id,
					source => $f->source(),
					genome => $genome,
					locus => $f->display_name,
					name => join('|',@notes),
					start => $f->start,
					end => $f->end,
					feature_id => $f->id,
					match => $value,
					matchtype => $header 
				};
				return $row;
			}
		}
	}
}

sub quick($){
	my $term = shift;
	my @results;
	my @termlist = split(/\s+/, $term);
	foreach $term (@termlist){
		my $query = "select fid from fattribute_to_feature where fattribute_value like '\%$term\%'";
		my $fids = query_get_fids($query);
		foreach my $fid (@$fids){
			$f = $db->get_feature_by_id($fid);
			push(@results, results_factory($f,$term,"Match"));
		}
	}
	my @sorted = sort {$b->{overlap} <=> $a->{overlap}} @results;
	return \@sorted;
}

###############Nemesys(Pubmed)###################

sub pubmed($){
	my $query = shift;
	my @results;
	my @features = $db->features(-types=>["CDS:Nemesys"]);
	foreach my $f (@features)
	{
		my @tagvalues = $f->get_tag_values('pubmedID');
		foreach $id (@tagvalues)
		{
			if($id =~ /$query/i)
			{
				my @notes = $f->get_tag_values('Location');
				my $overlap=scalar( $f->overlapping_features() );
				my $row = {
					overlap => $overlap,
					contig => $f->seq_id,
					source => $f->source(),
					genome => $genome,
					locus => $f->display_name,
					name => join('|',@notes),
					start => $f->start,
					end => $f->end,
					feature_id => $f->id,
					match => $id,
					matchtype => 'Pubmed ID'
				};
				push(@results,$row);
				last;
			}
		}
	}
	my @sorted = sort {$b->{overlap} <=> $a->{overlap}} @results;
	return \@sorted;
}


####################Kegg###################

sub keggorthology($)
{
	my $query = shift;
	my @results;
	my @features = $db->features(-types=>["CDS:KEGG"]);
	foreach my $f (@features)
	{
		my @tagvalues = $f->get_tag_values('dBxref');
		foreach $id (@tagvalues)
		{
			if($id =~ /$query/i)
			{
				my @notes = $f->get_tag_values('Note');
				my $overlap=scalar( $f->overlapping_features() );
				my $row = {
					overlap => $overlap,
					contig => $f->seq_id,
					source => $f->source(),
					genome => $genome,
					locus => $f->display_name,
					name => join('|',@notes),
					start => $f->start,
					end => $f->end,
					feature_id => $f->id,
					match => $id,
					matchtype => 'Kegg ID'
				};
				push(@results,$row);
				last;
			}
		}	
	}	

	my @sorted = sort {$b->{overlap} <=> $a->{overlap}} @results;
	return \@sorted;
}




################# InterProScan##################

sub interpro($)
{

	my $query = shift;
	my @results;


	my @features = $db->features(-types=>[
		"CDS:BlastProDom",
		"CDS:FPrintScan",
		"CDS:Gene3D",
		"CDS:HAMAP",
		"CDS:HMMPfam",
		"CDS:HMMPIR",
		"CDS:HMMSmart",
		"CDS:HMMTigr",
		"CDS:PatternScan",
		"CDS:ProfileScan",
		"CDS:SIGNAL-P_HMM",
		"CDS:SIGNAL-P_NN",
		"CDS:superfamily",
		"CDS:TMHMM"
		]);

	foreach my $f (@features)
	{
		my @tagvalues = $f->get_tag_values('DBxref');
		foreach $id (@tagvalues)
		{
			if ($id =~ /$query/i)
			{
				my @notes = $f->get_tag_values('Note');
				my $overlap=scalar( $f->overlapping_features() );
				my $row = {
					overlap => $overlap,
					contig => $f->seq_id,
					source => $f->source(),
					genome => $genome,
					locus => $f->display_name,
					name => join('|',@notes),	
					start => $f->start,
					end => $f->end,
					feature_id => $f->id,
					match => $id,
					matchtype => 'IprScan ID'
				};
				push(@results,$row);
				last;
			}	
		}
	}	
	if(! scalar @results){
		push(@results,{matchtype=>"Interpro ID",match=>"$query:no matches found"});
	}

	my @sorted = sort {$b->{overlap} <=> $a->{overlap}} @results;
	return \@sorted;
}


#################Uniprot/SwissProt############

sub swissprot($)
{
	my $query = shift;
	my @results;
	my @features = $db->features(-types=>["CDS:Uniprot/Swissprot"]);
	foreach my $f (@features){
		my @tagvalues = $f->get_tag_values('dBxref');
		foreach $id (@tagvalues){
			if ($id =~ /$query/i){
				my $overlap=scalar( $f->overlapping_features() );
				my $row = {
					overlap => $overlap,
					contig => $f->seq_id,
					source => $f->source(),
					genome => $genome,
					locus => $f->display_name,
					name => join('|',@tagvalues),
					start => $f->start,
					end => $f->end,
					feature_id => $f->id,
					match => $id,
					matchtype => 'Swissprot ID'
				};
				push(@results,$row);
				last;
			}
		}
	}	

	my @sorted = sort {$b->{overlap} <=> $a->{overlap}} @results;
	return \@sorted;
}
sub keggpathway($){
	my $query = shift;
	my @results;
	
	my @features = $db->features(-types=>["CDS:KEGG"]);
	
	foreach my $f (@features)
	{
		my @tagvalues = $f->get_tag_values('pathway');
	
		foreach $id (@tagvalues)
		{
			if ($id =~ /$query/i)
			{
				my $overlap=scalar( $f->overlapping_features() );
				my $row = {
					overlap => $overlap,
					contig => $f->seq_id,
					source => $f->source(),
					genome => $genome,
					locus => $f->display_name,
					name => join('|',@tagvalues),	
					start => $f->start,
					end => $f->end,
					feature_id => $f->id,
					match => $id,
					matchtype => 'KEGG Pathway'
				};
				push(@results,$row);
				last;
			}	
		}
	}
	my @sorted = sort {$b->{overlap} <=> $a->{overlap}} @results;
	return \@sorted;
}

sub protein($)
{
	my $query = shift;
	my @results;
	
	
	my @features = $db->features(-types=>["gene"]);
	foreach my $f (@features)
	{
		my @tagvalues = $f->get_tag_values('label');
		foreach $id (@tagvalues)
		{
			if ($id =~ /$query/i)
			{
				my $overlap=scalar( $f->overlapping_features() );
				my $row = {
					overlap => $overlap,
					contig => $f->seq_id,
					source => $f->source(),
					genome => $genome,
					locus => $f->display_name,
					name => join('|',@tagvalues),	
					start => $f->start,
					end => $f->end,
					feature_id => $f->id,
					match => $id,
					matchtype => 'Protein'
				};
				push(@results,$row);
				last;
			}	
		}
	}
	my @sorted = sort {$b->{overlap} <=> $a->{overlap}} @results;
	return \@sorted;
}


sub gene($$$$){
	my $keyword = shift;
	my $min = shift;
	my $max = shift;
	my $feature_type = shift;
	
	my @results;
	my @features;
	
	if($min!=0 || $max!=0)
	{	
		my $keyword = "select fid from fdata where ";
		if($min){$keyword .= "fstop - fstart >= $min";}
		if($min && $max){$keyword .= " and ";}
		if($max){$keyword .= "fstop - fstart <= $max";}
		my $fids = query_get_fids($keyword);
	
		foreach my $fid (@$fids)
		{
			my $ftr = $db->get_feature_by_id($fid);
			if($ftr->primary_tag() !~ /^$feature_type$/){ next;}
			push (@features,$ftr);
		}
	
	}
	
	else
	{
		@features = $db->features(-types=>[$feature_type]);
	}
	foreach my $f (@features)
	{
		my @tagvalues = $f->get_tag_values('label');
		foreach $id (@tagvalues)
		{
			if($id =~ /$keyword/i)
			{
				$id =~ s/_/ /g;
				my $overlap=scalar( $f->overlapping_features() );
				my $row = {
					overlap => $overlap,
					contig => $f->seq_id,
					source => $f->source(),
					genome => $genome,
					locus => $f->display_name,
					name => join('|',@tagvalues),
					start => $f->start,
					end => $f->end,
					feature_id => $f->id,
					match => $id,
					matchtype => 'gene'
				};
				push(@results,$row);
				last;
			}
		}
	}
	my @sorted = sort {$b->{overlap} <=> $a->{overlap}} @results;
	return \@sorted;
}

sub details($$$){ ## TO DO: make a template file for the html
	my ($contig,$Name,$f) = @_;
	$contig = $f->seq_id; # contig
	my $output .= '<div class="subfeature">';
	my @features = $db->overlapping_features(-refseq => $contig, -start => $f->start, -stop => $f->end );
	$Name=$f->display_name;
	$output .= "<table><tbody><tr><th colspan=3>Associated Features</th></tr><tr>";
	my $counter=0;
	foreach my $subftr (@features){
		if ($subftr->display_name =~ /^$contig$/){next;} # do not print contig
		if ($subftr->primary_id() == $f->primary_id()){next;}#same feature
		$output .= '<td><a href="/gb2/gbrowse_details/nemo?feature_id=' .
			$subftr->primary_id() . '">' . $subftr->display_name . '</a><ul>';
		my @tags = $subftr->get_all_tags();
		$output .= "<li><b>type:</b> " . join(', ',$subftr->type()) . "</li>";
		foreach my $tag (@tags){
			my $taginfo = join(', ',$subftr->get_tag_values($tag));
			if($tag =~ /^dBxref$/){
				$taginfo=parse_dbxref($taginfo);
			}
			if($tag =~ /^pathway$/){
			#	print STDERR "dbxref " . join('',$subftr->get_tag_values('dBxref'));
				$taginfo=parse_kegg_pathways($taginfo,join('',$subftr->get_tag_values('dBxref')));
			}
			$taginfo =~ s/_/ /g;
			$output .= "<li><b>$tag:</b> " . display_string($taginfo) . "</li>";
		
		}
		$output .= "</ul></td>";
		$counter = (1+$counter) % 4;
		if($counter>2){$output .= "</tr><tr>";}
	}
	$output .= "</tr></tbody></table>";
	$output .="</div><br/>";
	$output .= display_string(join('',$f->get_tag_values('label'))); 
	my $gourl = "http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term=";
	$output =~ s/(GO:\d+)/<a href="$gourl$1">$1<\/a>/g;
	$output;
}
sub parse_kegg_pathways($$){
	my ($pathway,$keggid) = @_;
	$keggid =~ s/^[^:]+://;
	print STDERR "Kegg ID $keggid";
	my $ecurl = "http://www.genome.jp/dbget-bin/www_bget?ec:";
	my $pathurl = "http://www.genome.jp/kegg-bin/show_pathway?ko";
	$pathway =~ s/\[EC:([^\]]+)\]/<a href="$ecurl$1">[EC:$1]<\/a>/g;
	$pathway =~ s/(\d{5})/<a href="$pathurl$1+$keggid">$1<\/a>/g;
	$pathway =~ s/\|/<br\/>/g;
	$pathway;
}

sub parse_dbxref($$$){
	my ($dbxref) = @_;
	my $iprurl = "http://www.ebi.ac.uk/interpro/IEntry?ac=IPR";
	$dbxref =~ s/(InterPro_ID:IPR)(\d+)/<a href="$iprurl$2">$1$2<\/a>/g;
	$dbxref;
}
sub display_string($){
	my $input = shift;
	#$input =~ s/_/ /g;
	$input =~ s/\|/, /g;
	return $input;
}
sub getcontigs ($){
	my $genome = shift;	
	my @values = ();
	my $db = gff_connect($genome);
	my @features = $db->features(-type=>"contig");
	foreach my $feature (@features){
		push(@values,sprintf("%s",$feature->name()));
	}
	return \@values;
}
sub gff_connect($){
	my $dbname = shift;
	return Bio::DB::GFF->new( -adaptor=>'dbi::mysqlopt',-dsn=>"dbi:mysql:$dbname",-user=>'gbrowse2',-pass=>'Brows3r2') or die "ERROR in CGPBase\n";
}
