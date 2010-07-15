#!/usr/bin/perl -w
use warnings;
use strict;
use CGI qw/:standard/;
use CGI qw/:cgi-lib/;
use lib '/var/www/nbase2010/lib'; #Perl Template Toolkit
use CGPBase;
use Template;
my $tpl = Template->new({INCLUDE_PATH => '/var/www/nbase2010/html/view',}) or  die "template error\n";
my $data = '';
print header('text/html');
my $vars = Vars;
$vars={%$vars, organisms=>CGPBase::getorganisms()};
my $page;
if (defined $$vars{'page'}) {
	$page=$$vars{'page'};
}
else {
	$page="about";
}
$$vars{current}='home';
if (exists($$vars{'raw'})){
	$tpl->process("$page.tt",$vars) or die $tpl->error();
}
else{
	$tpl->process("$page.tt",$vars,\$data) or die $tpl->error();
	$$vars{'content_for_layout'} = $data;
	$tpl->process('page.tt', $vars);
}
