#!/usr/bin/env perl
use strict;
use warnings;
my %cogid2gene;
my %gene2cogid;
open FH, "<$ARGV[0]" or die $!;
my @lines=<FH>;
close FH;
chomp @lines;
my $cogid="NONE";
foreach my $line (@lines){
	$line =~ s/^\s*(.*)\s*$/$1/;
	if($line=~/^$/){next;}
	if($line=~/^[_]+$/){$cogid="NONE";next;}
	my @words=split(/\s+/,$line);
	if($words[0]=~/^\[[^]]+\]/){
		$cogid=$words[1];
		@{$cogid2gene{$cogid}}=();
		next;
	}
	while($words[0]=~/^.*:$/){shift @words;}#organism code
	push @{$cogid2gene{$cogid}},@words;
	foreach my $gene (@words){
		if(!exists $gene2cogid{$gene}){@{$gene2cogid{$gene}}=();}
		push @{$gene2cogid{$gene}},$cogid;
	}
}
#foreach my $key (keys %cogid2gene){
#	printf "\t%s: %s\n",$key,join(', ',@{$cogid2gene{$key}});
#}
foreach my $key (keys %gene2cogid){
	printf "%s\t%s\n",$key,join(', ',@{$gene2cogid{$key}});
}
