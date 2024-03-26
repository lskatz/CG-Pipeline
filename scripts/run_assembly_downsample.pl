#!/usr/bin/env perl 

# run_assembly_downsample.pl: Downsample a set of reads
# Author: Lee Katz <lkatz@cdc.gov>
use strict;
use warnings;
use File::Basename qw/basename/;

print STDERR basename($0).": Running run_assembly_removeDuplicateReads.pl --nobin\n\n";
system("run_assembly_removeDuplicateReads.pl --nobin @ARGV");
die if $?;

