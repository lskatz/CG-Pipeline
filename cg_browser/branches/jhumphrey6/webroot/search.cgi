#!/usr/bin/perl -w
#use warnings;
use strict;
no strict "refs";
#use Bio::DB::GFF;
#use Data::Dumper;
use CGI qw/:standard/;
use CGI qw/:cgi-lib/;
use FindBin;
use lib ("$FindBin::Bin/../lib"); #Perl Template Toolkit
use CGPBase;
use Template;
use Data::Dumper;
my $tpl = Template->new({INCLUDE_PATH => "$FindBin::Bin/view",}) or  die "template error\n";
my $vars = Vars; # input variables from HTTP request
	
sub getpage($$){
	my ($template,$myvars) = @_;
	my $output = "";
	$tpl->process($template,$myvars,\$output) or die $tpl->error();
	return $output;
}
sub getfullpage($$){
	my ($template,$myvars) = @_;
	my $output = "";
	my $finalout = "";
	$tpl->process($template,$myvars,\$output) or die $tpl->error();
	$$myvars{'content_for_layout'} = $output;
	$tpl->process('page.tt', $myvars,\$finalout);
	return $finalout;
}
sub debug{
	print STDERR Dumper(@_);
}


#main#
$vars={%$vars,
	organisms=>CGPBase::getorganisms(),
	featuretypes=>CGPBase::gettypes(),
	sources=>CGPBase::getsources(),
	current=>'search'
};
print header('text/html');
if (defined $$vars{'sub'}) {
	my $action = $$vars{'sub'};
	if ($action =~ /^filters$/){
		print getpage("search/filters.tt",$vars);
	}
	if ($action =~ /^organism$/){
		my $contigs=();
		my $genome="unknown";
		if(defined $$vars{'genome'}){
			$genome=$$vars{'genome'};
			#$contigs=getlist($genome);
			$contigs=CGPBase::getcontigs($genome);
		}
		$vars={%$vars,genome=>$genome,contigs=>$contigs};
		print getpage("search/organism.tt",$vars);
	}
	elsif (grep (/^$action$/, qw/quick gene protein sig_peptide tmhmm/)){
		my $cgpcall = "CGPBase::$action";
		if($$vars{execute}){
			my @results = &$cgpcall($vars);
			@$vars{'results'}=@results;
			$$vars{'external_link'}=$ENV{REQUEST_URI} . "&external";
			if(defined $$vars{external}){print getfullpage('search/searchresults.tt', $vars);}
			else{print getpage("search/searchresults.tt",$vars);}
		}
		else{
			print getpage("search/$action.tt",$vars);
		}
	}
	elsif ($action eq 'blast'){
		if($$vars{execute}){
			if(!$$vars{'f_organisms'}){
				my $orgs=join('|',@{$$vars{'organisms'}});
				$$vars{'f_organisms'}=$orgs;
			}
			$$vars{'blastroot'}="$FindBin::RealBin/../etc/blast";
			my @results=CGPBase::blast($vars);
			@$vars{'results'}=@results;
			print getpage("search/blastresults.tt",$vars);
		}
		else{
			print getpage("search/blast.tt",$vars);
		}
	}
	elsif ($action eq 'organism_data'){
		$$vars{'organism_data'}=CGPBase::organism_data;
		$$vars{'current'}='organism_data';
		print getfullpage('organism_data.tt', $vars);
	}
	else {
	#	print STDERR "controller not found\n";
	}
}
else{#send the basic search forms
	print getfullpage('search/search.tt', $vars);
}
