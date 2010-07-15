#!/usr/bin/perl -w
use warnings;
use strict;
use CGI qw/:standard/;
use CGI qw/:cgi-lib/;
use Data::Dumper;
use lib '/var/www/nbase2010/lib'; #Perl Template Toolkit
use Template;
my $tpl = Template->new({INCLUDE_PATH => '/var/www/nbase2010/html/view',}) or  die "template error\n";
my $data = '';
print header('text/html');
my $vars = {current=>'download'};
my $dir;
opendir($dir,"download") or die "$!\n";
my @files=sort readdir($dir);
for my $file (@files){
	if(! -d $file){
		push(@{$$vars{'files'}},$file);
	}
}
$tpl->process("download.tt",$vars,\$data) or die $tpl->error();
$$vars{'content_for_layout'} = $data;
$tpl->process('page.tt', $vars);
