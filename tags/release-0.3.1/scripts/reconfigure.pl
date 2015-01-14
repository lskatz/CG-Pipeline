#!/usr/bin/env perl

# configure.pl: modify conf/cgpipelinerc
# Author: Jay Humphrey (jhumphrey6@gatech.edu)

package PipelineRunner;

use FindBin;
use File::Basename;
use Date::Format;
use File::Copy;

my $cfgfile="$FindBin::RealBin/../conf/cgpipelinerc";
if($ARGV[0]){$cfgfile=$ARGV[0];}
open(FH,"<$cfgfile") or die "Cannot open $cfgfile:$!\n";
my @content=<FH>;
close(FH);
print "$cfgfile:" . $#content . "\n";
my %config;
my $INSTALL_PATH=dirname dirname $FindBin::RealBin;
my @index;
foreach (@content){
	chomp;
	if(/^\s*([^=]+)\s*=\s*([^#]+)$/){
		my ($key,$val)=($1,$2);
		$key=~s/^\s*(\S*)\s*$/$1/;
		$val=~s/^\s*(\S*)\s*$/$1/;
		$config{$key} = $val;
		push(@index,$key);
	}
}
print "CG-Pipeline configuration\nCurrent settings appear in [brackets] - Press <Enter> to accept the current setting\n";
foreach my $key (@index){
	if(grep(/^$key$/,qw/pipeline_version pipeline_reference/)){next;}
	my $val=$config{$key};
	print "$key [$val] >:";
	my $input=<STDIN>;
	chomp $input;
	if($input!~/^$/){
		$config{$key}=$input;
	}
}
print "\nWrite new configuration? [Y/n]:";
my $cmd=<STDIN>; chomp $cmd; $cmd = uc $cmd;
if($cmd eq 'N'){print "Aborting\n";exit 0;}
my $save_old=$cfgfile . "." . time2str("%Y%m%d%H%M",time);
print "\nSaving old $cfgfile as $save_old\n";
copy $cfgfile, $save_old or die "Cannot write $save_old:$!\n";
print "\nWriting new configuration\n";
open(FH,">$cfgfile") or die "Cannot write $cfgfile:$!\n";
foreach my $key (@index){
	my $val=$config{$key};
	print FH "$key = $val\n";
}
close(FH);
print "\nDone\n";
