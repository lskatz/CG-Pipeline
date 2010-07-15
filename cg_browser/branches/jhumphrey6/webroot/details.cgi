#!/usr/bin/perl -w
use warnings;
use strict;
no strict "refs";
use Bio::DB::GFF;
use Data::Dumper;
use CGI qw/:standard/;
use CGI qw/:cgi-lib/;
use lib '/var/www/nbase2010/lib'; #Perl Template Toolkit
use CGPBase;
use Template;
my $tpl = Template->new({INCLUDE_PATH => '/var/www/nbase2010/html/view',}) or  die "template error\n";
my $data = {};
my $vars = Vars; # input variables from HTTP request


sub getfullpage($$){
	my ($template,$myvars) = @_;
	my $output = "";
	my $finalout = "";
	$tpl->process($template,$myvars,\$output) or die $tpl->error();
	$$myvars{'content_for_layout'} = $output;
	$tpl->process('details/page.tt', $myvars,\$finalout);
	return $finalout;
}

if($$vars{'file'}){
	if ($$vars{'file'} eq 'fasta'){
		print header('text/plain');
		CGPBase::fasta($$vars{'source'},$$vars{'fid'});
	}
}
else{
	print header('text/html');
	$data = CGPBase::details_page($$vars{'source'},$$vars{'fid'});
	print getfullpage("details/details.tt",$data);
}
