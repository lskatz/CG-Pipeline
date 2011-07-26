#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;
use File::Basename;

if(!defined $ARGV[0]){exit;}
my $makefile=File::Spec->rel2abs(shift);
while($makefile =~ /\.\./){$makefile =~ s/(\/[^\/]+\/\.\.)//;}
my $destdir=dirname($makefile);

open(IF,"<$makefile") or die "Cannot read $makefile:$!\n";
my @content=<IF>;
close(IF);
open(OF,">$makefile") or die "Cannot write $makefile:$!\n";
foreach (@content){
	s/^(\s*DESTDIR\s*=\s*)[^\s]+$/$1"$destdir"/;
	print OF;
}
close(OF);
