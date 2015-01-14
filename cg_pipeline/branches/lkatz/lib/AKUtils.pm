#!/usr/bin/env perl

# AKUtils.pm: utility subroutine library
# Author: Andrey Kislyuk (kislyuk@gatech.edu)

package AKUtils;
require 5.005;
my ($VERSION) = ('$Id$' =~ /,v\s+(\d+\S+)/o);

use strict;
use IPC::Open2;
use Cwd;
use List::Util qw(min max sum reduce shuffle);
use File::Basename;
use File::Spec;
use File::Copy;
use File::Temp ('tempdir');
use Sys::Hostname;
use Storable ('dclone'); # for deep copying
#use Bit::Vector;
use Data::Dumper; # useful for printing deep structures, e.g. to print %hash: print Data::Dumper->new([\%hash],[qw(hash)])->Indent(3)->Quotekeys(0)->Dump;
#use Fatal qw(:void open);
# Useful modules to look at:
# FindBin DynaLoader IO::Poll IO::Select IO::Socket Math::* Net::*
# algorithm::diff archive::tar authen::pam compress::* datemanip crypt::* event::* soap-lite svg term::* BerkeleyDB
# Perl tail recursion: @_=(arg1, arg2, arg3); goto &subroutine_name;

# TODO: allpairs method: all adjacent pairs in an array, also implement with checking for in-place modification
# TODO: seqwiz -> AKGenomeUtils (how to maintain processivity?)

use threads;
use threads::shared;
use Thread::Queue;

use Exporter;
our @ISA = "Exporter";
our @methods = qw(argmax compactArray concatFiles curSub filesDiffer flatten fullPathToExec nanmean safeGlob indexMfa getSeq getSubSeq printSeqsToFile lruQueue readPSL getFreeMem totalSize isFasta readMfa allIndexes startStopPos nanmin nanmax numEltsAtDepth alnum permuteRanges loadCSV intersectLength gaussian logmsg mktempdir loadConfig is_fastqPE);
our %EXPORT_TAGS = (all => [@methods]);
Exporter::export_ok_tags('all');

local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };

our %codon_table = (
	'TTT'=>'F', 'TTC'=>'F', 'TTA'=>'L', 'TTG'=>'L',
	'CTT'=>'L', 'CTC'=>'L', 'CTA'=>'L', 'CTG'=>'L',
	'ATT'=>'I', 'ATC'=>'I', 'ATA'=>'I', 'ATG'=>'M',
	'GTT'=>'V', 'GTC'=>'V', 'GTA'=>'V', 'GTG'=>'V',
	'TCT'=>'S', 'TCC'=>'S', 'TCA'=>'S', 'TCG'=>'S',
	'CCT'=>'P', 'CCC'=>'P', 'CCA'=>'P', 'CCG'=>'P',
	'ACT'=>'T', 'ACC'=>'T', 'ACA'=>'T', 'ACG'=>'T',
	'GCT'=>'A', 'GCC'=>'A', 'GCA'=>'A', 'GCG'=>'A',
	'TAT'=>'Y', 'TAC'=>'Y', 'TAA'=>'*', 'TAG'=>'*',
	'CAT'=>'H', 'CAC'=>'H', 'CAA'=>'Q', 'CAG'=>'Q',
	'AAT'=>'N', 'AAC'=>'N', 'AAA'=>'K', 'AAG'=>'K',
	'GAT'=>'D', 'GAC'=>'D', 'GAA'=>'E', 'GAG'=>'E',
	'TGT'=>'C', 'TGC'=>'C', 'TGA'=>'*', 'TGG'=>'W',
	'CGT'=>'R', 'CGC'=>'R', 'CGA'=>'R', 'CGG'=>'R',
	'AGT'=>'S', 'AGC'=>'S', 'AGA'=>'R', 'AGG'=>'R',
	'GGT'=>'G', 'GGC'=>'G', 'GGA'=>'G', 'GGG'=>'G',
);

our %seqs;
our (%seq_cache, %cdb_handle_cache);
our (@seq_cache_queue); # need a linkedhashmap type of thing
our ($max_seq_cache_items, $max_seq_cache_size) = (1024, 200000);

our ($cur_cdb_file, $cdb_pid, $cur_getseq_file, $seq_cache_size, $totsize_warned, $devel_size_loaded, $string_approx_loaded, $quiet, $leftover_gaussian);
local (*CDB_RH, *CDB_WH, *GETSEQ_RH);

END {
	close GETSEQ_RH if defined fileno GETSEQ_RH;
#	close CDB_WH if defined fileno CDB_WH;
#	close CDB_RH if defined fileno CDB_RH;
#	waitpid($cdb_pid, 0) if defined $cdb_pid;
}

sub logmsg {my $FH = $AKUtils::LOG || *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

sub mktempdir(;$) {
	my ($settings) = @_;
	my $tempdir_path = File::Spec->join(File::Spec->tmpdir(), (split("::",(caller(1))[3]))[1].".$$.XXXXX");
	my $tempdir = tempdir($tempdir_path, CLEANUP => !($$settings{keep}));
	return $tempdir;
}

sub commonSetup($) {
	# Accept a settings hash with options for setup
	# Return a settings hash with tmpdir, etc

	# Load common modules, include current path in INC and PATH, set up tmpdir, export useful symbols, parse options, etc.
}

sub dna2aa($) {
	my ($dna_seq) = @_;
	my $aminoacid_seq;
	for (my $i=0; $i<length($dna_seq); $i += 3) {
		my $codon = substr($dna_seq, $i, 3);
		$aminoacid_seq .= $codon_table{$codon};
		$aminoacid_seq .= "X" unless $codon_table{$codon};
	}
	return $aminoacid_seq;
}

# Loads data from file using Data::Dumper
sub loadData($) {
	my ($file) = @_;
	open(IN, '<', $file) or die("Unable to open file $file for reading: $!");
	my $ref;
	{	no strict 'vars';
		local $/;
		my $d = <IN>;
		$ref = eval $d;
		use strict 'vars';
	}
	close IN;
	return $ref;
}

# Saves data to a file using Data::Dumper
sub saveData($$) {
	my ($ref, $file) = @_;
	require Data::Dumper;
}

# Assumes the keys are unique. TODO: keyless load
# Expects header line, keys on column 0 unless otherwise specified (by column index or name)
# TODO: sanity checks
sub loadCSV($;$) {
	my ($file, $key) = @_;
	my %table; my $key_col;
	open(FH, '<', $file) or die("Could not open file $file for reading: ".$!);
	my $header = <FH>; chomp $header;
	my @fieldnames = split(/,/, $header);
	my %f; for (@fieldnames) { die("Empty or duplicate header field") if $f{$_} or $_ eq ''; $f{$_} = 1; }
	$key = 0 if not defined $key;
	for (0..$#fieldnames) { $key_col = $_ if $key eq $fieldnames[$_]; }
	$key_col = $key if $key =~ /^\d+$/;
	while (<FH>) {
		chomp;
		my @fields = split /,/;
		for (0..$#fieldnames) {
			$table{$fields[$key_col]}->{$fieldnames[$_]} = $fields[$_];
		}
	}
	close FH;
	return(\%table, \@fieldnames);
}

# TODO: optimize this
# TODO: setting: unbiased vs. biased estimator
# Can use Statistics::Descriptive.
sub stdev {
	my @list = @_;
	return undef if @list < 1;
#	$list = $$list[0] if @$list == 1 and ref($$list[0]) == 'ARRAY';
	return undef if @list < 1;
	my $mean = sum(@list) / @list;
#	return sqrt(sum(map { ($_-$mean)**2 } @$list) / (@$list-1));
	return sqrt(sum(map { ($_-$mean)**2 } @list) / @list);
}

# TODO: setting: averaging allowed
# Can use Statistics::Descriptive.
sub median {
	my @list = @_;
	return undef if @list < 1;
#	$list = $$list[0] if @$list == 1 and ref($$list[0]) == 'ARRAY';

	my @slist = sort {$a <=> $b} @list;
	my $count = scalar @list;
	if ($count % 2) {
		return $slist[int($count/2)];
	} else {
		return ($slist[$count/2] + $slist[$count/2 - 1]) / 2;
	}
}

# Return a value in the list such that $percentile percent of values in the list
# are below that value.
# TODO: setting: accept a user-defined sorting function and return index not value
sub percentile($$) {
	my ($list, $percentile) = @_;
	die("argument percentile undefined") if not defined $percentile;
	die("illegal value of argument percentile") if $percentile<0 or $percentile>100;
	return undef if not defined $list or @$list < 1;

	my @slist = sort {$a <=> $b} @$list;
	my $count = scalar @$list;
	return $slist[min(int($count*$percentile/100), $#slist)];
}

# can also use PDL
# Box-Muller transform (http://www.taygeta.com/random/gaussian.html)
sub gaussian($$) {
	my ($mean, $variance) = @_;
	if (defined $leftover_gaussian) {
		my $y = $leftover_gaussian;
		undef $leftover_gaussian;
		return $y;
	}
	my ($x1, $x2, $w, $y1, $y2);
	$variance = 1 unless defined $variance;

	do {
		$x1 = 2.0 * rand() - 1.0;
		$x2 = 2.0 * rand() - 1.0;
		$w = $x1 * $x1 + $x2 * $x2;
	} while ($w >= 1.0);

	$w = sqrt((-2.0 * log($w)) / $w);
	$y1 = ($x1 * $w) * $variance + $mean;
	$y2 = ($x2 * $w) * $variance + $mean;
	$leftover_gaussian = $y2;
	return $y1;
}

sub argmax {
	my ($hash) = @_;
	my ($best, $key, $val);

	my $bestval = -1e999;
	while (($key,$val) = each(%$hash)) {
		if ($val > $bestval) {
			$bestval = $val;
			$best = $key;
		}
	}
	return $best;
}

# http://www.perlmonks.org/?node_id=371228
# Usage:
#	my $iter = combinations( 3 => ['a' .. 'f'] );
#	while ( my @c = $iter->() ) { print "@c\n"; }
sub combinations($$) {
	my ($num, $arr) = @_;
	no strict 'refs';
	return sub { return }
		if $num == 0 or $num > @$arr;
	use strict 'refs';
	my @pick;

	return sub {
		return @$arr[ @pick = ( 0 .. $num - 1 ) ]
			unless @pick;

		my $i = $#pick;
		$i-- until $i < 0 or $pick[$i]++ < @$arr - $num + $i;
		return if $i < 0;

		@pick[$i .. $#pick] = $pick[$i] .. $#$arr;

		return @$arr[@pick];
	};
}

# TODO: optimize me
# Usage: @list_of_lists = permutations(@list)
sub permutations(@) {
	my (@set) = @_;
	return [@set] if @set < 2;
	my @perms;
	foreach my $i (0..$#set) {
		my @tail_perms = permutations(@set[0..$i-1, $i+1..$#set]);
		foreach my $subperm (@tail_perms) {
			push(@perms, [$set[$i], @$subperm]);
		}
	}
	return @perms;
}

# Input: { a => [b, c], d => [e, f], ... }
# Output: [ { a => b, d => e }, { a => b, d => f }, ... ]
sub permuteRanges($) {
	my ($ranges) = @_;
	my $permutations = [{}];
	foreach my $attrib (keys %$ranges) {
		my $new_permutations;
		foreach my $value (@{$$ranges{$attrib}}) {
			my $new_set = dclone($permutations); # deep copy
			foreach my $p (@$new_set) {
				$$p{$attrib} = $value;
			}
			push(@$new_permutations, @$new_set);
		}
		$permutations = $new_permutations;
	}
	return $permutations;
}

sub sampleWithoutReplacement($$) {
	my ($set, $sample_size) = @_;
	die("Sample size should be a positive integer") if $sample_size < 1 or int($sample_size) != $sample_size;
	die("Sample size should not exceed set size") if $sample_size > @$set;
	my @sample;
	my %taken_indices;
	for (1..$sample_size) {
		my $index;
		do { $index = int(rand(scalar(@$set))) } while defined $taken_indices{$index};
		$taken_indices{$index} = 1;
		push(@sample, $$set[$index]);
	}
	return wantarray ? @sample : \@sample;
}

# Remove null, zero elements from array
# Input, output: array ref
sub compactArray {
	my $array = (@_ == 1 ? $_[0] : \@_);
	my @new_array;
	foreach my $element (@$array) { push(@new_array, $element) if $element; }
	return wantarray ? @new_array : \@new_array;
}

# Remove zero elements from array
# Input, output: array ref
sub compactArrayWithZeros {
	my $array = (@_ == 1 ? $_[0] : \@_);
	my @new_array;
	foreach my $element (@$array) { push(@new_array, $element) if defined $element; }
	return wantarray ? @new_array : \@new_array;
}

# TODO: make this use File::Temp
sub concatFiles($;$) {
	my ($files, $tmp_dir_or_file) = @_;
	die("No files supplied") unless @$files > 0;
	$tmp_dir_or_file = ($ENV{TMPDIR} or "/tmp") unless $tmp_dir_or_file;
	if (-d $tmp_dir_or_file) {
		my $rand_name;
		do { $rand_name = substr(rand(), 2, 4) } while -e "$tmp_dir_or_file/$$.$rand_name";
		$tmp_dir_or_file .= "/$$.$rand_name";
	}
	open(OUT, '>', $tmp_dir_or_file) or die("Could not open file $tmp_dir_or_file for writing: ".$!);
	foreach my $file (@$files) {
		open(FH, '<', $file) or die("Could not open file $file for reading: ".$!);
		print OUT while <FH>;
		close FH;
	}
	close OUT;
	return $tmp_dir_or_file;
}

sub curSub() {
	return (caller(1))[3];
}

sub filesDiffer($$) {
# TODO: use Algorithm::Diff
# @f1=<FH>; $diff = Algorithm::Diff->new(\@f1, \@f2); next   if  $diff->Same();
	my ($file1, $file2) = @_;
	open(FH1, '<', $file1) or die("Could not open file $file1 for reading: ".$!);
	open(FH2, '<', $file2) or die("Could not open file $file2 for reading: ".$!);
	while (1) {
		my $line1 = <FH1>;
		my $line2 = <FH2>;
		last unless defined $line1 and defined $line2;
		return 1 if $line1 ne $line2;
	}
	return 1 if <FH1> or <FH2>;
	return 0;
}

# Traverse array of arrays down to scalars and make it a flat array of those scalars.
# optimize me
sub flatten($) {
	my ($array) = @_;
	my @flattened;

	for (ref($array) eq 'ARRAY' ? @$array : ref($array) eq 'HASH' ? values %$array : $array) {
		push(@flattened, (ref($_) ? flatten($_) : $_));
	}
	return @flattened;
}

# If argument is an executable in the current path, returns the full path to it, otherwise returns undef
# arguments: warn_on_error
sub fullPathToExec($;$) {
	my ($executable,$settings) = @_;
	my $fullpath;
	for ("", split(/:/, $ENV{PATH})) {
    my $path=$_."/".$executable;
		if (-x $path && -f $path) { $fullpath = File::Spec->rel2abs($path); last; }
	}
  if(! -x $fullpath){
	  my $errStr="Error finding full path to executable ($executable)";
    logmsg $errStr if($$settings{warn_on_error});
    die $errStr if(!$$settings{warn_on_error});
  }
	return $fullpath;
}

sub safeGlob($) {
	my ($regexp) = @_;
	my $dir;
	if ($regexp =~ /\//) { ($dir, $regexp) = ($regexp =~ /\A(.*\/)([^\/]+)\Z/); }
	$regexp .= "\$" unless $regexp =~ /\$$/;
	$regexp = "^".$regexp unless $regexp =~ /^\^/;
	my (@files);
	local (*DIR);

	$dir ||= ".";
	$regexp ||= ".*";
	opendir(DIR, $dir) or return;
	@files = grep { /$regexp/ } readdir(DIR);
	closedir(DIR);
	foreach my $file (@files) {
		(undef($file), next) if $file eq "." or $file eq "..";
		$file = $dir.$file if $dir;
	}
	return @{compactArray(\@files)};
}

sub md5sum($) {
	my ($filename) = @_;
	my $md5sum;
	die("$filename: file not found") unless -f $filename;
	if (`which md5sum 2>/dev/null`) {
		$md5sum = `md5sum $filename`; $md5sum =~ /^([^\s]+)/; $md5sum = $1;
	} else {
		require Digest;
		open(FH, '<', $filename) or die("Could not open file $filename for reading: ".$!);
		my $md5obj = Digest->new('MD5');
		$md5obj->addfile(\*FH);
		$md5sum = $md5obj->hexdigest();
		close FH;
	}
	return $md5sum;
}

sub cdb_indexMfa($;$$) {
	require Digest;
	my ($filename, $seqs, $quiet) = @_;
	$seqs = \%seqs if @_ < 2;
	my ($seq_name, $start_fpos, $end_fpos, $seq_len, $last_line, $no_cached_index, $md5sum, $old_md5sum);

	print curSub().": Checksumming $filename... \n" unless $quiet;
	open(FH, '<', $filename) or die("Could not open file $filename for reading: ".$!);
	my $md5obj = Digest->new('MD5');
	$md5obj->addfile(\*FH);
	$md5sum = $md5obj->hexdigest();
	close FH;

	if (-f $filename.'.md5') {
		open(FH, '<', $filename.'.md5') or die("Could not open file $filename.md5 for reading: ".$!);
		$old_md5sum = <FH>;
		chomp $old_md5sum;
		close FH;
	}

	print curSub().": Loading stored index for $filename... \n" if $md5sum eq $old_md5sum and not $quiet;
	if ($md5sum ne $old_md5sum) {
		warn(curSub().": Warning: Index checksum mismatch in $filename. Re-indexing.\n") if defined $old_md5sum and not $quiet;
		print curSub().": Indexing $filename... \n" unless $quiet;
		my $invoke_string = "cdbfasta $filename"; $invoke_string .= ">/dev/null 2>&1" if $quiet;
		system($invoke_string); die("Error $? executing \"$invoke_string\". Stopped") if $?;
		open(FH, '>', $filename.'.md5') or die("Could not open file $filename.md5 for writing: ".$!);
		print FH $md5sum."\n";
		close FH;
	}

	open(IN, "cdbyank -l $filename.cidx|") or die("Could not run \"cdbyank -l $filename.cidx\": ".$!);
	while (<IN>) {
		/^([^\t]+)\t(\d+)$/ or die; # sequence name, length
		$$seqs{$1} = [$filename, $2];
	}
	close IN;
	return $seqs;
}

# Load MFA index for $filename, append it to %$seqs, return the md5sum for the MFA file from which the index was created
=head1
sub doIndexMfaFile($$) {
	my ($filename, $seqs) = @_;
	my ($seq_name, $start_fpos, $end_fpos, $seq_len);

	print curSub().": loading $filename\n";
	open(FH, "< ".$filename) or die(curSub().": Could not open file $filename for reading: ".$!);
	while (<FH>) {
		if (/^\>[\s]*([^\s]+)/) {
			my $new_seq_name = $1;
			if (defined $seq_name) {
				warn(curSub().": duplicate sequence name \"$seq_name\" in $filename") if defined $$seqs{$seq_name};
				die(curSub().": empty sequence \"$seq_name\" in $filename") unless $start_fpos < $end_fpos;
				$$seqs{$seq_name} = [$filename, $start_fpos, $end_fpos, $seq_len];
			}
			$start_fpos = tell FH;
			$seq_name = $new_seq_name; $seq_len = 0;
		} else {
			chomp; s/\s//g; #s/[^ATGCNatgcn]//g;
			$seq_len += length;
		}
		$end_fpos = tell FH;
	}
	warn(curSub().": duplicate sequence name \"$seq_name\" in $filename") if defined $$seqs{$seq_name};
	die(curSub().": empty sequence \"$seq_name\" in $filename") unless $start_fpos < $end_fpos;
	$$seqs{$seq_name} = [$filename, $start_fpos, $end_fpos, $seq_len];
	seek FH, 0, 0; # back to start
	my $md5obj = Digest->new("MD5");
	$md5obj->addfile(\*FH);
	my $md5sum = $md5obj->hexdigest();
	close FH;

#	print curSub().": done indexing $filename\n";
	if (open(FH, '>', $filename.'.faindex')) {
		foreach my $seq_name (keys %$seqs) {
			print FH "$seq_name\t${$$seqs{$seq_name}}[1]\t${$$seqs{$seq_name}}[2]\t${$$seqs{$seq_name}}[3]\n";
		}
		print FH $md5sum."\n";
		close FH;
		print curSub().": saved mfa index in $filename.faindex\n";
	} else {
		print STDERR curSub().": failed to save mfa index in $filename.faindex\n";
	}
}
=cut


# TODO: the only thing to optimize is index creation, saving and loading, pack it with pack() and gzip
# TODO: if write to cur dir fails, write to global cache dir
# TODO: corrupt or wrong faindex will pollute %seqs before being discarded

# Input: name of file with multi-fasta sequences, ref to hash to append to
# Output: hash(key->sequence name, value->[filename, seek position of first char in sequence, seek pos of last char, real length in bp])
sub indexMfa($;$$) {
	my ($filename, $seqs, $quiet) = @_;
	$seqs = \%seqs if @_ < 2;

	print curSub().": Checksumming $filename... \n" unless $quiet;
	my $md5sum = md5sum($filename);

	my $saved_md5sum;
	if (-f $filename.".faindex") {
		open(FH, '<', $filename.".faindex") or warn(curSub().": Could not open file $filename.faindex for reading: ".$!);
		while (<FH>) {
			/^([^\t]+)\t([\d]+)\t([\d]+)\t([\d]+)$/o or (/^([\w]+)$/, $saved_md5sum = $1, last) or return;
			$$seqs{$1} = [$filename, $2, $3, $4];
		}
		chomp $saved_md5sum;
	}

	if (-f $filename.".faindex" and $md5sum eq $saved_md5sum) {
		print curSub().": Loaded stored index for $filename\n" unless $quiet;
	} else {
		print curSub().": Indexing $filename... \n" unless $quiet;

		my ($seq_name, $start_fpos, $end_fpos, $seq_len);
		open(FH, "< ".$filename) or die("Could not open file $filename for reading: ".$!);
		while (<FH>) {
			if (/^\>[\s]*([^\s]+)/) {
				my $new_seq_name = $1;
				if (defined $seq_name) {
					warn(curSub().": duplicate sequence name \"$seq_name\" in $filename") if defined $$seqs{$seq_name};
					die("Empty sequence \"$seq_name\" in $filename") unless $start_fpos < $end_fpos;
					$$seqs{$seq_name} = [$filename, $start_fpos, $end_fpos, $seq_len];
				}
				$start_fpos = tell FH;
				$seq_name = $new_seq_name; $seq_len = 0;
			} else {
				chomp; s/\s//g; #s/[^ATGCNatgcn]//g;
				$seq_len += length;
			}
			$end_fpos = tell FH;
		}
		warn(curSub().": duplicate sequence name \"$seq_name\" in $filename") if defined $$seqs{$seq_name};
		die("Empty sequence \"$seq_name\" in $filename") unless $start_fpos < $end_fpos;
		$$seqs{$seq_name} = [$filename, $start_fpos, $end_fpos, $seq_len];
		seek FH, 0, 0; # back to start
		my $md5obj = Digest->new("MD5");
		$md5obj->addfile(\*FH);
		my $md5sum = $md5obj->hexdigest();
		close FH;

		if (open(FH, '>', $filename.'.faindex')) {
			foreach my $seq_name (keys %$seqs) {
				print FH "$seq_name\t${$$seqs{$seq_name}}[1]\t${$$seqs{$seq_name}}[2]\t${$$seqs{$seq_name}}[3]\n";
			}
			print FH $md5sum."\n";
			close FH;
			print curSub().": saved mfa index in $filename.faindex\n";
		} else {
			warn curSub().": failed to save mfa index in $filename.faindex\n";
		}
	}
	return $seqs;
}

sub getSeq($;$) {
	my ($seq_name, $seqs) = @_;
	die("Error: sequence name not supplied") unless $seq_name;
	$seqs ||= \%seqs;
#print "getseq: retrieving $seq_name from ".join(":",keys %$seqs)."\n";

	my $retr_seq;
	my ($filename, $start_fpos, $end_fpos) = @{$$seqs{$seq_name}};
	die("Error: cannot find sequence seek info. Stopped") unless $filename and defined $start_fpos and $end_fpos;

#	my $csize; $csize = sum(map { $csize += length $$_ } values %seq_cache);

	my $cache_line = $seq_cache{$filename.':'.$seq_name};
	return $cache_line if defined $cache_line;

# EVICTION
# estimate new sequence to be ($end_fpos - $start_fpos) bytes
# totalSize is expensive so estimate
	$seq_cache_size += ($end_fpos - $start_fpos);

	my $shrink_seq_cache = 0;
	my $cur_seq_not_added = 1;
	my $must_free_mem_for_seq = max($end_fpos - $start_fpos - getFreeMem(), 0);
# must execute once to add tag initially!

#while ($cur_seq_not_added or $seq_cache_size > max($max_seq_cache_size,
# getFreeMem())) { - WRONG - must be "size of seq to add exceeds free mem"
# evict items from queue until it's small enough...
	my $line_to_evict = lruQueue(\@seq_cache_queue, $max_seq_cache_items - $shrink_seq_cache, $filename.':'.$seq_name);
	if ($line_to_evict) {
# TODO: check that the length below equals $end_fpos - $start_fpos
		$seq_cache_size -= length(${$seq_cache{$line_to_evict}});
		undef ${$seq_cache{$line_to_evict}};
		delete $seq_cache{$line_to_evict};
	}
	die("Internal error") unless $cur_seq_not_added or $line_to_evict;
	$cur_seq_not_added = 0;
	$shrink_seq_cache++;
#}
# END EVICTION


#print "getseq: not in cache\n";
	# NB: this must be line buffered, but Perl does that by default so we're fine
	if ($cur_getseq_file ne $filename) {
#print "getseq: reopening FH for $filename\n";
		close GETSEQ_RH if defined fileno GETSEQ_RH;
		open(GETSEQ_RH, '<', $filename) or die("Could not open file $filename for reading: ".$!);
		$cur_getseq_file = $filename;
	}
#print "getseq: $seq_name in $filename rlen ".($end_fpos - $start_fpos)."\n";

	seek(GETSEQ_RH, $start_fpos, 0);
	read(GETSEQ_RH, $retr_seq, $end_fpos - $start_fpos);
	$retr_seq =~ s/[\s\n]//g;

	$seq_cache{"$filename:$seq_name"} = \$retr_seq;
#print "getseq: done retrieving $seq_name\n";
	return \$retr_seq;
}

# TODO: make this cache aware
sub getSubSeq($$$;$) {
	my ($seq_name, $start, $length, $seqs) = @_;
	my $seq = getSeq($seq_name, $seqs);
	return substr($$seq, $start, $length);
}

sub cdb_printSeqsToFile($$;$) {
	my ($seqs_to_print, $to_file, $seqs) = @_;
	die("Error: sequence name not supplied") unless @$seqs_to_print > 0;
	$seqs ||= \%seqs;
	my %seqs_by_filename;

	foreach my $seq_name (@$seqs_to_print) {
		my ($filename, $seq_len) = @{$$seqs{$seq_name}};
		die("Error: cannot find sequence seek info. Stopped") unless defined $filename and defined $seq_name;
		push(@{$seqs_by_filename{$filename}}, $seq_name);
	}

	foreach my $filename (keys %seqs_by_filename) {
		# TODO: make processive together with getseq
		open(OUT, "|cdbyank $filename.cidx >> $to_file");
		print OUT join("\n", @{$seqs_by_filename{$filename}});
		close OUT;
	}
#	system("cdbyank $filename.cidx -a '$seq_name' >> $to_file");
}

# Output: reference to a string
# TODO: implement cdb handle cache (is it necessary though, with seq cache?)
# TODO: this is insanely slow
sub cdb_getSeq($;$) {
	my ($seq_name, $seqs) = @_;
	die("Error: sequence name not supplied") unless $seq_name;
	$seqs ||= \%seqs;

#print "getseq: retrieving $seq_name\n";
	my $retr_seq;
	my ($filename, $seq_len) = @{$$seqs{$seq_name}};
	die("Error: cannot find sequence seek info. Stopped") unless defined $filename and defined $seq_name;

	my $cache_line = $seq_cache{$filename.':'.$seq_name};
	my $line_to_evict = lruQueue(\@seq_cache_queue, $max_seq_cache_items, $filename.':'.$seq_name);
	die("Internal error") if $cache_line and $line_to_evict; # if it was in the cache, there should be nothing to evict
	delete $seq_cache{$line_to_evict} if $line_to_evict;
	return $cache_line if defined $cache_line;
#print "getseq: not in cache\n";
	#open(IN, "cdbyank $filename.cidx -a '$seq_name'|"); while(<IN>) { next if /^\>/; chomp; $retr_seq .= $_; } close IN;
	if ($cur_cdb_file ne $filename) {
#print "getseq: setting up handle\n";
		close CDB_WH if defined fileno CDB_WH;
		close CDB_RH if defined fileno CDB_RH;
		waitpid($cdb_pid, 0) if defined $cdb_pid;
		$cdb_pid = open2(\*CDB_RH, \*CDB_WH, "cdbyank $filename.cidx"); # or die...
		$cur_cdb_file = $filename;
	}
#print "getseq: requesting $seq_name\n";
	# NB: this must be line buffered, but Perl does that by default so we're fine
	print CDB_WH $seq_name."\n";

	while (<CDB_RH>) {
		next if /^\>/; # TODO: sanity check if debug: see if seq_name is present in fasta header
		last if /^\n/;
		chomp;
		$retr_seq .= $_; # TODO: more sanity checks, avoid deadlock in case of error
	}

	$seq_cache{"$filename:$seq_name"} = \$retr_seq;
	return \$retr_seq;
}

# TODO: this needs more work
sub isFasta($) {
	my ($file) = @_;
	my $char;
	open(FH, "< $file") or die("Could not open file $file for reading: ".$!);
	read(FH, $char, 1); close FH;
	if ($char eq ">") { return 1; } else { return 0; }
}

# Input: queue, max items in it, label of item to add to queue
# Output: label of item to be deleted, if any
# TODO: variable strategies, dynamic queue #items, optimize with a DLL
sub lruQueue($$$) {
	my ($queue, $max_queue_items, $item_to_add, $strategy) = @_;
	die("Bad arguments") unless $queue and $max_queue_items and defined $item_to_add;
	my $item_exists_at;
	foreach my $i (0..$#$queue) {
		if ($$queue[$i] eq $item_to_add) { $item_exists_at = $i; last; }
	}
	if (defined $item_exists_at) {
		splice(@$queue, $item_exists_at, 1);
		unshift(@$queue, $item_to_add);
#my @new_queue = ($item_to_add);push(@new_queue, @$queue[0..$item_exists_at-1]) if $item_exists_at > 0;push(@new_queue, @$queue[$item_exists_at+1..$#$queue]) if $item_exists_at < $#$queue;@$queue = @new_queue;
	} else {
		unshift(@$queue, $item_to_add);
	}

	if (@$queue > $max_queue_items) {
		die("Internal error") if @$queue > $max_queue_items+1;
		return pop(@$queue);
	} else {
		return undef;
	}
}

=head1
   1. matches - Number of bases that match that aren't repeats
   2. misMatches - Number of bases that don't match
   3. repMatches - Number of bases that match but are part of repeats
   4. nCount - Number of 'N' bases
   5. qNumInsert - Number of inserts in query
   6. qBaseInsert - Number of bases inserted in query
   7. tNumInsert - Number of inserts in target
   8. tBaseInsert - Number of bases inserted in target
   9. strand - '+' or '-' for query strand. For translated alignments, second '+'or '-' is for genomic strand
  10. qName - Query sequence name
  11. qSize - Query sequence size
  12. qStart - Alignment start position in query
  13. qEnd - Alignment end position in query
  14. tName - Target sequence name
  15. tSize - Target sequence size
  16. tStart - Alignment start position in target
  17. tEnd - Alignment end position in target
  18. blockCount - Number of blocks in the alignment (a block contains no gaps)
  19. blockSizes - Comma-separated list of sizes of each block
  20. qStarts - Comma-separated list of starting positions of each block in query
  21. tStarts - Comma-separated list of starting positions of each block in target
=cut
# $clusters = {key: db_seq, value: [[align1, ... alignN], [align1, ... alignN], ...]
sub readPSL($;$) {
	my ($filename, $clusters) = @_;
	$clusters ||= {};
	open(FH, '<', $filename) or die("Could not open file $filename for reading: ".$!);
	$_ = <FH>;
	warn(curSub().": File $filename doesn't look like a PSL file") unless /^psLayout version 3/;
	while (<FH>) {
		/^(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\+|-)\t([^\t]+)\t(\d+)\t(\d+)\t(\d+)\t([^\t]+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t([\d\,]+)\t([\d\,]+)\t([\d\,]+)$/
			or next;
		my ($n_match, $n_mism, $repm, $n_count, $q_nins, $q_bpins, $t_nins, $t_bpins, $strand, $q_name, $q_size, $q_start, $q_end, $t_name, $t_size, $t_start, $t_end, $block_cnt, $block_sizes, $q_starts, $t_starts)
			= ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21);

		push(@{$$clusters{$t_name}}, []);
		my $cur_cluster = ${$$clusters{$t_name}}[$#{$$clusters{$t_name}}];

		my ($i1, $i2, $i3);
		for (1..$block_cnt) {
			my ($next_i1, $next_i2, $next_i3) = (index($block_sizes, ',', $i1), index($q_starts, ',', $i2), index($t_starts, ',', $i3));
			my ($cur_blocksize, $cur_qstart, $cur_tstart)
				= (substr($block_sizes, $i1, $next_i1-$i1), substr($q_starts, $i2, $next_i2-$i2), substr($t_starts, $i3, $next_i3-$i3));
			($i1, $i2, $i3) = ($next_i1+1, $next_i2+1, $next_i3+1);

push(@$cur_cluster, [$q_name, $q_start, $q_end]);
		}
	}
	close FH;
	return $clusters;
}

# Usage: @indexes = allIndexes(string, substring)
# or @fuzzy_indexes = allIndexes(string, substring, #insertions, #deletions, #substitutions)
# WARNING: string::approx can produce wrong results in many situations, especially for short strings
# TODO: case insensitive, warn on min length/complexity w/amatch
# TODO: matches start at only certain position ranges or modulos
sub allIndexes($$;$$$) {
	my ($string, $substring, $in, $del, $sub) = @_;
	my @matches;
	if ($in or $del or $sub) {
		unless ($string_approx_loaded) {
			foreach my $dir (@INC) {
				if (-f "$dir/String/Approx.pm") {
					require String::Approx;
					$string_approx_loaded = 1;
				}
			}
			die("Unable to load String::Approx") unless $string_approx_loaded;
		}
		# string::approx is weird and can get stuck, limit the number of matches
		my ($nummatches, $maxmatches) = (0, 1000);
		my $pm = -length($substring);
		while (($pm = String::Approx::aindex($substring, ['I'.$in.'D'.$del.'S'.$sub.' initial_position='.($pm+length($substring))], $string)) != -1) {
			$nummatches++;
			(warn("Warning: Possible runaway match engine, indexes truncated to $maxmatches"), last) if $nummatches > $maxmatches;
			push(@matches, $pm);
		}
	} else {
		my $pm = -1;
		push(@matches, $pm) while ($pm = index($string, $substring, $pm+1)) != -1;
	}
	return wantarray ? @matches : \@matches;
}

# returns total size of structure and all its refs, in bytes, or undef if Devel::Size not installed
sub totalSize($) {
	my ($ptr) = @_;
	return undef unless $ptr;
	return Devel::Size::total_size($ptr) if $devel_size_loaded;
	foreach my $dir (@INC) {
		if (-f "$dir/Devel/Size.pm") {
			require Devel::Size;
			$devel_size_loaded = 1;
			return Devel::Size::total_size($ptr);
		}
	}
	warn(curSub().": Warning: Devel::Size not found. Unable to determine memory consumption\n") unless $totsize_warned or $quiet;
	$totsize_warned = 1;
	return undef;
}

# on Linux, returns estimated amount of free memory (MemFree + Cached), in bytes, or undef if unknown
sub getFreeMem() {
	my $free_mem;
	open(FH, '<', '/proc/meminfo') or return undef;
	while (<FH>) {
		if (/^MemFree.+?(\d+)\s+kB/) {
			$free_mem = $1;
		} elsif (/^Cached.+?(\d+)\s+kB/) {
			$free_mem += $1; last;
		}
	}
	close FH;
	return ($free_mem * 1024 or undef);
}

sub nanmean {
	my $array = ((@_ == 1 and ref($_[0]) eq 'ARRAY') ? $_[0] : \@_);
	my $size = (defined $$array[0] ? 1 : 0);
	my $sum = reduce { $size++ if defined $b; $a + $b } @$array;
	return undef if $size == 0;
	return $sum / $size;
}

sub nanmin {
	my $array = ((@_ == 1 and ref($_[0]) eq 'ARRAY') ? $_[0] : \@_);
	return undef unless $array;
	shift @$array while @$array and not defined $$array[0];
	my $min = reduce { defined $b ? $a < $b ? $a : $b : $a } @$array;
	return $min;
}

sub nanmax {
	my $array = ((@_ == 1 and ref($_[0]) eq 'ARRAY') ? $_[0] : \@_);
	return undef unless $array;
	shift @$array while @$array and not defined $$array[0];
	my $max = reduce { defined $b ? $a > $b ? $a : $b : $a } @$array;
	return $max;
}

# Compute the number of items/pairs at given traverse depth of a tree of arrays or hashes, e.g. numEltsAtDepth([[1,2],[3,4]], 1) = 4
sub numEltsAtDepth($$) {
	my ($ref, $depth) = @_;
	my @subrefs = ($ref);
	my $total;
	@subrefs = map { ref($_) eq 'ARRAY' ? @$_ : ref($_) eq 'HASH' ? values %$_ : die "Scan depth exceeds structure depth" } @subrefs while $depth--;
	map { $total += (ref($_) eq 'ARRAY' ? @$_ : ref($_) eq 'HASH' ? values %$_ : die "Scan depth exceeds structure depth") } @subrefs;
	return $total;
}

# TODO: error checking
# Usage: $seqs = readMfa($seqfile); foreach $seqname (keys %$seqs) { $seq = $$seqs{$seqname}; }
# Usage: (values %{readMfa($seqfile)})[0]
# TODO: return second array with ordering
sub readMfa($;$) {
	my ($mfa_file, $settings) = @_;
	open(FH, '<', $mfa_file) or die("Could not open file $mfa_file for reading: ".$!);
	my %seqs;
	my ($cur_seq_hdr, $seq);
	while (<FH>) {
		if (/^\>\s*(.+)/) {
			$seqs{$cur_seq_hdr} = $seq if $cur_seq_hdr;
			undef $seq;
			$cur_seq_hdr = ($$settings{first_word_only} ? (split /\s/, $1)[0] : $1);
		} else {
			chomp; 
      s/\s//g unless($$settings{keep_whitespace});
      $seq .= uc $_;
		}
	}
	close FH;
	$seqs{$cur_seq_hdr} = $seq if $cur_seq_hdr;
	die("Error: No sequences found in $mfa_file") unless %seqs;
	return \%seqs;
}

sub readGff($) {
	my ($gff_file) = @_;
	my @features;
	open(GFF, '<', $gff_file) or die "Error opening file $gff_file for reading: $!";
	while (<GFF>) {
		next if /^\#/;
		chomp;
		my @l = split /\t/;
		my %feature;
		foreach my $field (qw(seqid source type start end score strand phase attributes)) {
			my $value = shift @l;
			$feature{$field} = $value if $value ne '.';
		}
		# NB: Though the standard says attributes are in the form "attr1=value1;attr2=value2",
		# many files contain "attr1 value1;attr2 value2"
		foreach my $attrib (split /;/, $feature{attributes}) {
			$attrib =~ /^(.+)=(.+)$/ or next;
			foreach my $tag (qw(ID Name Alias Parent Target Gap Derives_from Note Dbxref Ontology_term)) {
				# NB: these can be separated by commas
				$feature{$1} = $2 if $1 eq $tag;
			}
		}
		push(@features, \%feature);
	}
	return \@features;
}

# TODO: this is a low-level parser, break the ugly out of readPSL and make that a high-level parser
sub readPSL2($) {
	my ($psl_file) = @_;
	my %hits;
	open(PSL, '<', $psl_file) or die "Error opening file $psl_file for reading: $!";
	<PSL> for 1..5; # skip header
	while (<PSL>) {
		chomp;
		my @l = split /\t/;
		my %hit;
		foreach my $field (qw(matches misMatches repMatches nCount qNumInsert qBaseInsert tNumInsert tBaseInsert strand qName qSize qStart qEnd tName tSize tStart tEnd blockCount blockSizes qStarts tStarts)) {
			my $value = shift @l;
			$hit{$field} = $value;
		}
		push(@{$hits{$hit{tName}}}, \%hit);
	}
	close PSL;
	return \%hits;
}

# TODO: GC% dependence
# Sequence: Single sequence to mutate
# Indels: Target # indels to introduce
# Substs (TODO)
# Rate: Indel rate per Kchar (either rate or indels must be defined)
# Subregions: A list of subregions in the form [{start=>, end=>}, {start=>, end=>}, ...] to insert errors into
# TODO: record coordinates must be adjusted when producing indels upstream of them
sub mutateSequence($$$;$$$) {
	my ($seq, $indels, $substs, $rate, $subregions, $whiteout) = @_;
	my @record;
	die("Cannot accept both rate and fixed target") if $rate and ($indels or $substs);
	die("invalid rate") if defined $rate and ($rate >= 1 or $rate <= 0);
	die("invalid target") if defined $indels and (int($indels) != $indels or $indels < 0);
	die("invalid target") if defined $substs and (int($substs) != $substs or $substs < 0);
	die("invalid subregions") if defined $subregions and ref($subregions) ne 'ARRAY';
	die("invalid whiteout") if defined $whiteout and ($whiteout < 1 or $whiteout > 1000);

	$rate = 1/10000 unless $rate or $indels or $substs;

	if ($rate and not ($indels or $substs)) {
		$indels = int(length($seq) * $rate);
		die("Internal error computing mutations from rate") if $indels > 10000;
		return ($seq, \@record) if $indels < 1;
	}

	$subregions ||= [[1, length($seq)-1]];

	MUT: while (1) {
		my $pos = int(rand(length($seq)));
		CDS: foreach my $region (@$subregions) {
			my ($start, $end) = ($$region{lo}, $$region{hi});
			die("Error reading subregions") unless $start and $end;
#			die("Subregion coordinate $end exceeds sequence length ".(length($seq)+) if $end > length($seq); # this check needs adjustment
			if ($pos > $start+$whiteout and $pos < $end-$whiteout) {
				#print "position $pos admitted in REGION $start $end\n";
				if (rand() < 0.5) { # insert
					substr($seq, $pos, 1) .= substr("ATGC", int(rand(length("ATGC"))), 1);
					push(@record, [$pos, 'I']);
				} else { # delete
					substr($seq, $pos, 1) = "";
					push(@record, [$pos, 'D']);
				}
				$indels--; last MUT if $indels == 0;
				last CDS;
			}
		}
	}
	return ($seq, \@record);
}

sub alnum {
	my ($i);
	my ($len1, $len2) = (length($a), length($b));
	for ($i = 0; ($i < $len1) && ($i < $len2); ++$i) {
		my $c1 = substr($a, $i, 1);
		my $c2 = substr($b, $i, 1);
		($c1 =~ /^\d/o) || ($c2 =~ /^\d/o) || ($c1 ne $c2) and last;
	}
	my $a_r = ($i < $len1) ? substr($a, $i) : "";
	my $b_r = ($i < $len2) ? substr($b, $i) : "";
	my ($a_n, $a_s) = ($a_r =~ /^(\d+)(.*)$/o);
	my ($b_n, $b_s) = ($b_r =~ /^(\d+)(.*)$/o);
	return (defined($a_n) && defined($b_n)) ?
	(($a_n <=> $b_n) || ($a_s cmp $b_s)) : ($a cmp $b);
}

sub intersectLength($$$$) {
	my ($start1, $end1, $start2, $end2) = @_;
	my $ilen;
	die("Internal error") unless $start1 and $end1 and $start2 and $end2;
	die("Internal error") unless $start1 <= $end1 and $start2 <= $end2;
#	($start1, $end1) = ($end1, $start1) if $start1 > $end1;
#	($start2, $end2) = ($end2, $start2) if $start2 > $end2;
	if ($start1 < $start2) {
		if ($end1 < $end2) {
			$ilen = $end1 - $start2 + 1;
		} else {
			$ilen = $end2 - $start2 + 1;
		}
	} else {
		if ($end2 < $end1) {
			$ilen = $end2 - $start1 + 1;
		} else {
			$ilen = $end1 - $start1 + 1;
		}
	}
	$ilen = 0 if $ilen < 0;
	return $ilen;
}

# Returns tRNAscanSE predictions as {seqname => [{lo =>, hi =>, strand =>, trna_type =>, anticodon =>}, ...], ...}
# Settings: optional: tse_exec: executable name for tRNAscanSE, tempdir: path to temporary directory
sub gettRNAscanSEPredictions($;$) {
	my ($seqs, $settings) = @_;
	die("Internal error: no input supplied") if ref($seqs) ne 'HASH';
	$$settings{tempdir} ||= AKUtils::mktempdir($settings);

	# tRNAscanSE mangles sequence names, so save and restore them
	my $i; my (%renamed_seqs, %seqname_orig2working, %seqname_working2orig);
	foreach my $seqname (keys %$seqs) {
		$i++;
		$seqname_orig2working{$seqname} = 'seq'.$i;
		$seqname_working2orig{'seq'.$i} = $seqname;
		$renamed_seqs{$seqname_orig2working{$seqname}} = $$seqs{$seqname};
	}
	my $input_file_full = "$$settings{tempdir}/trnascan-se.in.fasta";
	AKUtils::printSeqsToFile(\%renamed_seqs, $input_file_full);

	$$settings{trnascan_modeltype}||='bacteria';
	$$settings{tse_exec}||='tRNAscan-SE';
	$$settings{tse_exec} = AKUtils::fullPathToExec($$settings{tse_exec},{warn_on_error=>1});
	if(!$$settings{tse_exec}){
		logmsg "WARNING: could not locate tRNAscan-SE and so I will not run it. In the future, you can set the option tse_exec to point to the full path of tRNAscan-SE.";
		return {};
	}

  my $tse_opts="";
	$tse_opts .= "-q" if($$settings{quiet});
	$tse_opts .= "" if $$settings{trnascan_modeltype} eq 'eukaryotic';
	$tse_opts .= " -G" if $$settings{trnascan_modeltype} eq 'general';
	$tse_opts .= " -O" if $$settings{trnascan_modeltype} eq 'organellar';
	$tse_opts .= " -B" if $$settings{trnascan_modeltype} eq 'bacteria';
	$tse_opts .= " -A" if $$settings{trnascan_modeltype} eq 'archaea';
  $tse_opts .= " -Q"; # overwrite, to make it more automatic

        my $command="$$settings{tse_exec} $tse_opts -o '$$settings{tempdir}/trnascan-se.out' '$input_file_full'";
	logmsg "Running tRNAscanSE on $input_file_full with command\n  $command";
	system($command);
	die("Error running $$settings{tse_exec} with command $command: $!") if $?;

	return loadtRNAscanSEPredictions("$$settings{tempdir}/trnascan-se.out", \%seqname_working2orig);
}

sub loadtRNAscanSEPredictions($;$) {
	my ($predfile, $seqname_remap) = @_;
	open(PRED, '<', $predfile) or die("Could not open file $predfile for reading: ".$!);
	my $header .= <PRED> for 1..3;

	my %predictions;
	while (<PRED>) {
		chomp;
		my @l = split /\t/;
		s/^\s+// for @l; s/\s+$// for @l;
		my ($seqname, $index, $begin, $end, $trna_type, $anticodon, $intron_begin, $intron_end, $cove_score) = @l;
		$begin = int($begin); $end = int($end);
		my $strand = ($begin < $end ? '+' : '-');
		$seqname = $$seqname_remap{$seqname} if defined $seqname_remap;
		my $codon = reverse($anticodon); $codon =~ tr/ATGC/TACG/;
		push(@{$predictions{$seqname}}, {type => 'tRNA',
            seqname => $seqname, lo => min($begin, $end), hi => max($begin, $end), strand => $strand,
			trna_type => $trna_type, trna_anticodon => $anticodon, trna_codon_recognized => $codon,
            start => $begin, stop => $end,
			intron_begin => $intron_begin, intron_end => $intron_end, cove_score => $cove_score,
			predictor => 'tRNAscanSE'});
	}
	close PRED;
	return \%predictions;
}

# TODO: this might still be wrong - compare with bitmap-based
# TODO: check partial (end-of-sequence) orfs
# TODO: add options to discard ORFs with excessive Ns, require start (Met), accept partial (end-of-seq) orfs
sub findOrfs($$) {
	my ($seqs, $min_len) = @_;
	my ($start_codons, $stop_codons) = ("ATG|GTG", "TAA|TAG|TGA"); # FIXME: use startCodons()
	my ($rs_start_codons, $rs_stop_codons) = (scalar(reverse($start_codons)), scalar(reverse($stop_codons)));
	$rs_start_codons =~ tr/ATGC/TACG/; $rs_stop_codons =~ tr/ATGC/TACG/;
	my %orfs;

	while (my ($seqname, $seq) = each(%$seqs)) {
		warn("Sequence $seqname contains spaces or newlines") if $seq =~ /(\n|\s)/;
		warn("Sequence $seqname contains characters other than ATGCN") unless $seq =~ /^[ATGCN]+$/;

		foreach my $i (0, 1, 2) { # frame
			pos($seq) = $i;
			while ($seq =~ /\G((?:.{3})+?(?:$stop_codons))/go) { # forward strand
				my $orf = $1;
	#			if ($orf =~ /($start_codons)((?:.{3})+)$/o and length($2)+3 >= $min_len) { # contains at least one in-frame start
				if (length($orf)+3 >= $min_len) {
					$orfs{$seqname}->{$i}->{pos($seq)} = pos($seq)-length($orf)+1;
				}
			}
		}

		my $rc_seq = reverse($seq); $rc_seq =~ tr/ATGC/TACG/;
		foreach my $i (0, 1, 2) { # frame
			pos($rc_seq) = ((length($rc_seq) % 3) - $i) % 3;
			while ($rc_seq =~ /\G((?:.{3})+?(?:$stop_codons))/go) {
				my $orf = $1;
	#			if ($orf =~ /($start_codons)((?:.{3})+)$/o and length($2)+3 >= $min_len) { # contains at least one in-frame start
				if (length($orf)+3 >= $min_len) {
					$orfs{$seqname}->{$i.'R'}->{length($rc_seq)-pos($rc_seq)+1} = length($rc_seq) - pos($rc_seq) + length($orf);
				}
			}
		}
	}
	return \%orfs;
}

=head1
Extract prokaryotic (uninterrupted) ORFs from a nucleotide sequence.
Output: {seq1name => {frame => {stop_coord => {start=>N, stop=>N, seq=>..., aa_seq=>...}}}}
Settings:
	min_len: minimum length of ORFs to report. Default is 90 nt.
	need_seq: output ORF nucleotide sequence in the "seq" field. Default is false.
	need_aa_seq: output ORF amino acid sequence in the "aa_seq" field (starting at the reported start coordinate). Default is false.
x	confirm_start: Require a Met codon to be in the ORF, and modify the meaning of min_len to mean
x		the distance from the most upstream Met codon to the end of the ORF.
	get_truncated_orfs: Report ORFs truncated by the start or end of the sequence. The confirm_start setting is
		ignored for these ORFs; the min_len setting is in effect unless min_partorf_len is set. The coordinate
		closest to the truncated end is reported as the first nucleotide of a full in-frame codon on that end
		(i.e. for sequence start, it's 1, 2, or 3). The fields truncated_lo or truncated_hi are set accordingly.
x	min_partorf_len: Require truncated ORFs to be at least this length. Default is equal to min_len.
	orftype: if set to start2stop, 
The output fields are:
	ustop: the first nucleotide downstream of the stop codon of the preceding in-frame ORF
	start: the most upstream nucleotide of the first Met codon in the ORF
	stop: the most downstream nucleotide of the stop codon in the ORF
	truncated_upstream: set to true if the ORF may be truncated upstream (no in-frame stop codon precedes the ORF)
	truncated_downstream: set to true if the ORF is truncated downstream (sequence terminates before the ORF's stop codon)
	seq: nucleotide sequence from start to stop
	aa_seq: nucleotide sequence from start to stop
=cut
sub findOrfs2($;$) {
	my ($seqs, $settings) = @_;
	if (ref($seqs) ne 'HASH') { $seqs = {'seq1' => $seqs}; }
	logmsg "Extracting ORFs from ".keys(%$seqs)." sequences...";
	my ($start_codons, $stop_codons) = ("ATG|GTG", "TAA|TAG|TGA"); # FIXME: use startCodons()
	my ($rs_start_codons, $rs_stop_codons) = (scalar(reverse($start_codons)), scalar(reverse($stop_codons)));
	$rs_start_codons =~ tr/ATGC/TACG/; $rs_stop_codons =~ tr/ATGC/TACG/;
	$$settings{min_len} = 90 unless defined $$settings{min_len};
	my %orfs;
	while (my ($seqname, $seq) = each(%$seqs)) {
		die("Sequence $seqname contains spaces or newlines") if $seq =~ /(\n|\s)/;
		warn("Sequence $seqname contains characters other than ATGCN") unless $seq =~ /^[ATGCN]+$/;

		foreach my $frame (0, 1, 2) {
			my $cur_orf_start = $frame;
			my $cur_pos;
			for ($cur_pos = $frame; $cur_pos < length($seq)-2; $cur_pos += 3) { # walk codon-by-codon
				my $codon = substr($seq, $cur_pos, 3);
				if ($codon =~ /^(?:$stop_codons)$/ or $cur_pos+3 >= length($seq)-2) { # codon is a stop codon or this is the end of sequence
					my $orf_len = $cur_pos - $cur_orf_start + 3;
					if ($orf_len < $$settings{min_len}) {
						$cur_orf_start = $cur_pos+3; next;
					} elsif (not($$settings{get_truncated_orfs}) and $cur_pos+3 >= length($seq)-2 and $codon !~ /^(?:$stop_codons)$/) {
						$cur_orf_start = $cur_pos+3; next;
					}
					
					my $truncated_upstream = 1 if $cur_orf_start == $frame;
					my $truncated_downstream = 1 if $cur_pos+3 >= length($seq)-2 and $codon !~ /^(?:$stop_codons)$/;
					my $orf_nt_seq = substr($seq, $cur_orf_start, $orf_len);
					
					my $start_offset = findStartInORF($orf_nt_seq, $start_codons);
					if ($start_offset < 0 or $orf_len-$start_offset < $$settings{min_len}) {
						$cur_orf_start = $cur_pos+3; next;
					}
					if ($$settings{orftype} eq "start2stop") {
						$orf_nt_seq = substr($seq, $cur_orf_start+$start_offset, $orf_len-$start_offset);
					}

					my $stop = $cur_pos+3;
					$orfs{$seqname}->{$frame}->{$stop} = {
						ustop => $cur_orf_start+1,
						start => $cur_orf_start+1+$start_offset,
						stop => $cur_pos+3,
						length => $orf_len,
						strand => '+'};
					$orfs{$seqname}->{$frame}->{$stop}->{truncated_upstream} = 1 if $truncated_upstream;
					$orfs{$seqname}->{$frame}->{$stop}->{truncated_downstream} = 1 if $truncated_downstream;
					$orfs{$seqname}->{$frame}->{$stop}->{seq} = $orf_nt_seq if $$settings{need_seq};
					$orfs{$seqname}->{$frame}->{$stop}->{aa_seq} = dna2aa($orf_nt_seq) if $$settings{need_aa_seq};
					$orfs{$seqname}->{$frame}->{$stop}->{lo} = $orfs{$seqname}->{$frame}->{$stop}->{start};
					$orfs{$seqname}->{$frame}->{$stop}->{hi} = $orfs{$seqname}->{$frame}->{$stop}->{stop};
					$cur_orf_start = $cur_pos+3;
				}
			}
		}
		
		foreach my $frame (0, 1, 2) { # Reverse frame
			my $cur_orf_start = length($seq) - ((length($seq)-$frame) % 3) - 1;
			my $cur_pos;
			for ($cur_pos = $cur_orf_start; $cur_pos >= 2; $cur_pos -= 3) { # walk codon-by-codon
				my $codon = substr($seq, $cur_pos-2, 3);
				if ($codon =~ /^(?:$rs_stop_codons)$/ or $cur_pos-3 < 0) { # codon is a stop codon
					my $orf_len = $cur_orf_start - $cur_pos + 3;
					if ($orf_len < $$settings{min_len}) {
						$cur_orf_start = $cur_pos-3; next;
					} elsif (not($$settings{get_truncated_orfs}) and $cur_pos-3 < 0 and $codon !~ /^(?:$rs_stop_codons)$/) {
						$cur_orf_start = $cur_pos-3; next;
					}

					my $truncated_upstream = 1 if $cur_orf_start == length($seq) - ((length($seq)-$frame) % 3) - 1;
					my $truncated_downstream = 1 if $cur_pos-3 < 0 and $codon !~ /^(?:$rs_stop_codons)$/;

					my $orf_nt_seq = substr($seq, $cur_pos-2, $orf_len);
					die("Internal error") unless length($orf_nt_seq) % 3 == 0;
					$orf_nt_seq = reverse($orf_nt_seq); $orf_nt_seq =~ tr/ATGC/TACG/;

					my $start_offset = findStartInORF($orf_nt_seq, $start_codons);
					if ($start_offset < 0 or $orf_len-$start_offset < $$settings{min_len}) {
						$cur_orf_start = $cur_pos-3; next;
					}
					if ($$settings{orftype} eq "start2stop") {
						$orf_nt_seq = substr($seq, $cur_pos-2, $orf_len-$start_offset);
						$orf_nt_seq = reverse($orf_nt_seq); $orf_nt_seq =~ tr/ATGC/TACG/;
					}

					die("Internal error") unless length($orf_nt_seq) % 3 == 0;
					my $stop = $cur_pos-1;
					$orfs{$seqname}->{$frame.'R'}->{$stop} = {
						ustop => $cur_orf_start+1,
						start => $cur_orf_start+1-$start_offset,
						stop => $cur_pos-1,
						length => $orf_len,
						strand => '-',
					};
					$orfs{$seqname}->{$frame.'R'}->{$stop}->{truncated_upstream} = 1 if $truncated_upstream;
					$orfs{$seqname}->{$frame.'R'}->{$stop}->{truncated_downstream} = 1 if $truncated_downstream;
					$orfs{$seqname}->{$frame.'R'}->{$stop}->{seq} = $orf_nt_seq if $$settings{need_seq};
					$orfs{$seqname}->{$frame.'R'}->{$stop}->{aa_seq} = dna2aa($orf_nt_seq) if $$settings{need_aa_seq};
					$orfs{$seqname}->{$frame.'R'}->{$stop}->{lo} = $orfs{$seqname}->{$frame.'R'}->{$stop}->{stop};
					$orfs{$seqname}->{$frame.'R'}->{$stop}->{hi} = $orfs{$seqname}->{$frame.'R'}->{$stop}->{start};
					$cur_orf_start = $cur_pos-3;
				}
			}
		}
	}
	return \%orfs;
}

# Internal for findOrfs
sub findStartInORF($$) {
	my ($seq, $start_codons) = @_;
	for (my $cur_pos = 0; $cur_pos < length($seq)-2; $cur_pos += 3) {
		return $cur_pos if substr($seq, $cur_pos, 3) =~ /^(?:$start_codons)$/;
	}
	return -1;
}

# Run and parse SignalP output.
# Both the runner and the parser in BioPerl were found inadequate.
# Usage: $predictions = getSignalPPredictions($seqs, {signalp_org_type => 'gram-', ...})
#     $seqs = {seq1name => seq1, ...}
#     $predictions = {seq1name => [ {type=>..., ...}, ...], ...}
sub getSignalPPredictions($$) {
	my ($seqs, $settings) = @_;
	die("No input sequences supplied") if ref($seqs) ne 'HASH';
	die("Setting signalp_org_type must be set to gram-, gram+ or euk") if $$settings{signalp_org_type} !~ /^(gram\-|gram\+|euk)$/;
	die("Invalid value of setting signalp_pred_method, must be one of nn (neural networks), hmm (hidden Markov models), nn+hmm")
		if defined $$settings{signalp_pred_method} and $$settings{signalp_pred_method} !~ /^(nn|hmm|nn\+hmm)$/;
	$$settings{signalp_trunc_length} = 70 if not defined $$settings{signalp_trunc_length};
  $$settings{numcpus}||=AKUtils::getNumCPUs();
  #$$settings{numcpus}=4; # DEBUG
	die("Invalid value of setting signalp_trunc_length, must be an integer between 0 and 1000")
		if int($$settings{signalp_trunc_length}) != $$settings{signalp_trunc_length}
	        or $$settings{signalp_trunc_length} < 0 or $$settings{signalp_trunc_length} > 1000;
	$$settings{signalp_exec} ||= AKUtils::fullPathToExec('signalp');
	die("SignalP executable not found") unless -x $$settings{signalp_exec};

	logmsg "Running SignalP on ".values(%$seqs)." sequences using $$settings{numcpus} threads...";

	$$settings{tempdir} ||= AKUtils::mktempdir($settings);

	# See man signalp for option documentation
	my $signalp_opts = "-t $$settings{signalp_org_type}";
	$signalp_opts .= " -f summary"; # TODO: support full and short formats
	# TODO: support graphics
	$signalp_opts .= " -method $$settings{signalp_pred_method}" if defined $$settings{signalp_pred_method};
	$signalp_opts .= " -trunc $$settings{signalp_trunc_length}";
	
  my $Q=Thread::Queue->new;
  my @thr;
  for(0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(\&signalPWorker,$signalp_opts,$Q,$settings);
  }
	my $predictions = {};
	my $i;
	foreach my $seqname (sort alnum keys %$seqs) {
		$i++;
    next if(!$$seqs{$seqname});
    $Q->enqueue({$seqname=>$$seqs{$seqname}});
		logmsg "Enqueued $i proteins" if $i % 100 == 0;
    #last if($i>60); # DEBUG
	} $|++; $|--; # flush output
  $Q->enqueue(undef) for(@thr);
  while((my $pending=$Q->pending)>$$settings{numcpus}){
    logmsg "Waiting on $pending seqs to be processed";
    sleep 10;
  }
  for(0..$$settings{numcpus}-1){
    my $t=$thr[$_];
    logmsg "Waiting on thread ".$t->tid;
    my $thrPred=$t->join;
    $predictions={%$predictions,%$thrPred};
  }
	logmsg "Processed $i proteins, done";
	return $predictions;
}

sub signalPWorker{
  my ($signalp_opts,$Q,$settings)=@_;
  my $tid=threads->tid;
  my $tempdir=$$settings{tempdir};
  my $predictions = {};
  my $i=0;
  while(defined(my $seq=$Q->dequeue)){
    $i++;
    my ($seqid,$sequence)=each(%$seq);
		if (length($sequence) < 20 or length($sequence) > 2000) {
			logmsg "Peptide sequence $seqid length ".length($sequence)." out of bounds (20-2000 residues)";
			next if(length($sequence) < 20);
      # this next part runs if it is in the >2000aa realm
      logmsg "  I truncated $seqid to 2000 residues just for SignalP";
      $sequence=substr($sequence,0,2000);
		}
    my $signalpInfile="$tempdir/signalp.TID$tid.$i.in.fasta";
    my $signalpOutfile="$tempdir/signalp.TID$tid.$i.out";
		printSeqsToFile({$seqid => $sequence}, $signalpInfile) or die "Could not print seqs to file $signalpInfile";
		my $invoke_string = "$$settings{signalp_exec} $signalp_opts < $signalpInfile > $signalpOutfile";
		system($invoke_string);
		die("Error running signalp: $!") if $?;

		my $result = loadSignalPPredictions($signalpOutfile);
		$predictions = {%$predictions, %$result};
  }
  return $predictions;
}

sub loadSignalPPredictions($) {
	my ($predfile) = @_;
	my %predictions;
	open(PRED, '<', $predfile) or die("Could not open file $predfile for reading: ".$!);
	my $pred_text;
	{ local $/; $pred_text = <PRED>; }
	close PRED;
	my @pred_blocks = split(/^\-+$/, $pred_text);

	foreach my $pred_block (@pred_blocks) {
		$pred_block =~ /^>(.+)$/m or die("Internal error while parsing $predfile");
		my $name = $1;
		if ($pred_block =~ /^SignalP-NN result:(.+?)^$/ms) {
			my $nn_pred_block = $1;
			while($nn_pred_block =~ /^\s*(max. C|max. Y|max. S|mean S|D)\s+(.+)$/gm) {
				my %prediction = (type => 'nn');
				$prediction{measure} = $1;
				my @l = split(/\s+/, $2);
				for (qw(position value cutoff decision)) {
					$prediction{$_} = shift @l;
				}
				($prediction{start}, $prediction{end}) = split(/-/, $prediction{position});
				push(@{$predictions{$name}}, \%prediction);
			}
		}
		if ($pred_block =~ /^SignalP-HMM result:(.+?)^$/ms) {
			my $hmm_pred_block = $1;
			my %prediction = (type => 'hmm');
			$hmm_pred_block =~ /^Prediction: (.+)$/m;
			$prediction{decision} = $1;
			$hmm_pred_block =~ /^Signal peptide probability: (.+)$/m;
			$prediction{probability} = $1;
			$hmm_pred_block =~ /^Max cleavage site probability:\s+([^\s]+)\s+between pos.\s+([^\s]+)\s+and\s+([^\s]+)\s*$/m;
			($prediction{CLP}, $prediction{start}, $prediction{end}) = ($1, $2, $3);
			push(@{$predictions{$name}}, \%prediction);
		}
	}
	return \%predictions;
}

sub getGlimmer3Predictions($;$) {
	my ($seqs, $settings) = @_;
	die("Internal error: no input supplied") if ref($seqs) ne 'HASH';

	logmsg "Preparing to run glimmer3 on ".scalar(keys %$seqs)." sequences...";

	$$settings{tempdir} ||= AKUtils::mktempdir($settings);

	# Save and restore sequence names
	my $i; my (%renamed_seqs, %seqname_orig2working, %seqname_working2orig);
	foreach my $seqname (keys %$seqs) {
		$seqname_orig2working{$seqname} = $seqname;
		$seqname_working2orig{$seqname} = $seqname;
		$renamed_seqs{$seqname_orig2working{$seqname}} = $$seqs{$seqname};
		$i++;
    next;
		$seqname_orig2working{$seqname} = 'seq'.$i;
		$seqname_working2orig{'seq'.$i} = $seqname;
		$renamed_seqs{$seqname_orig2working{$seqname}} = $$seqs{$seqname};
	}

  if(-e "$$settings{tempdir}/glimmer3.predict"){
    logmsg "Found Glimmer3 output file at $$settings{tempdir}/glimmer3.predict. Not going to run it again";
    return loadGlimmer3Predictions("$$settings{tempdir}/glimmer3.predict", $seqs, \%seqname_working2orig);
  }

	my $glimmer_infile = "$$settings{tempdir}/glimmer3_in.fasta";
	printSeqsToFile(\%renamed_seqs, $glimmer_infile);

	my $longorfs_infile = $glimmer_infile;
	if (scalar(keys %$seqs) > 1) {
		$longorfs_infile = "$$settings{tempdir}/glimmer3_longorfs_in.fasta";
		open(FH, '>', $longorfs_infile) or die "Could not open file $longorfs_infile for writing: ".$!;
		# TODO: CHECK - spacer might bias model
		my $seqs = join('TTAGTTAGTTAG', values %$seqs); # sequences separated by all-frame stop spacer
		$seqs =~ s/(.{80})/$1\n/g;
		$seqs .= "\n" unless $seqs =~ /\n$/;
		print FH ">glimmer3_longorfs_in\n$seqs";
		close FH;
	}

	my ($longorfs_opts, $glimmer_opts);
	# Sequences are assumed to be non-circular if there's more than one of them
	$longorfs_opts = "--linear" if $$settings{linear_genome} or scalar(keys %$seqs) > 1;
	$glimmer_opts = " --linear" if $$settings{linear_genome} or scalar(keys %$seqs) > 1;
	$glimmer_opts .= " --rbs_pwm $$settings{gl3_rbs_pwm_file}" if $$settings{gl3_rbs_pwm_file};
	$glimmer_opts .= " --gene_len $$settings{gl3_min_gene_len}" if $$settings{gl3_min_gene_len};
	$glimmer_opts .= " --threshold $$settings{gl3_calling_threshold}" if $$settings{gl3_calling_threshold};
	$glimmer_opts .= " --max_olap $$settings{gl3_max_overlap}" if $$settings{gl3_max_overlap};
	$glimmer_opts .= " --extend" unless $$settings{gl3_no_extended_predicts};

  my $exec=AKUtils::fullPathToExec("long-orfs");
	system("$exec $longorfs_opts --no_header --cutoff 1.15 '$longorfs_infile' '$longorfs_infile.longorfs' 2>'$$settings{tempdir}/glimmer3.log'");
	die("Error running long-orfs with command\n  long-orfs $longorfs_opts --no_header --cutoff 1.15 '$longorfs_infile' '$longorfs_infile.longorfs' 2>'$$settings{tempdir}/glimmer3.log'\nError was $!") if $?;
	system("extract -t '$longorfs_infile' '$longorfs_infile.longorfs' > '$$settings{tempdir}/glimmer3.train'");
	die("Error running extract: $!") if $?;
	system("build-icm -r '$$settings{tempdir}/glimmer3.icm' < '$$settings{tempdir}/glimmer3.train'");
	die("Error running build-icm: $!") if $?;
	logmsg "Running glimmer3 on $glimmer_infile...";
	system("glimmer3 $glimmer_opts '$glimmer_infile' '$$settings{tempdir}/glimmer3.icm'"
		. " '$$settings{tempdir}/glimmer3' 2>>'$$settings{tempdir}/glimmer3.log'");
	die("Error running glimmer3: $!") if $?;

	return loadGlimmer3Predictions("$$settings{tempdir}/glimmer3.predict", $seqs, \%seqname_working2orig);
}

=head1
NB: Handling of predictions which extend past the end of the sequence
In gmhmmp, a prediction extending past the end is denoted by "<1", "<2", "<3" or ">[end-1]" etc. positions. Coordinates start at 1.
In glimmer3, it's denoted by a coordinate extending by 1, 2 or 3 nt beyond the edge of the sequence. Coordinates start at 1.
The objective here is to preserve length %3 == 0, but the resulting coordinates differ.
To regularize this behavior, I add/subtract 3 to/from the glimmer3 coordinates, remove </> from the gmhmmp coordinates,
and set the extended_lo or extended_hi flags in the prediction hash.
With Glimmer, this also requires the $seqs hash, since we need to know the sequences' lengths.
=cut
sub loadGlimmer3Predictions($;$$) {
	my ($predfile, $seqs, $seqname_remap) = @_;
	open(PRED, '<', $predfile) or die("Could not open file $predfile for reading: ".$!);
	my %predictions;
	my $cur_seq = 'unknown';
	while (<PRED>) {
		if (/^\>\s*(.+?)\s*$/) {
			$cur_seq = $1;
			$cur_seq = $$seqname_remap{$cur_seq} if(defined $seqname_remap && defined($$seqname_remap{$cur_seq}));
		} else {
			my @line = split /\s+/;
			my ($gene_id, $start, $end, $frame, $raw_score) = @line;
			my $strand = ($frame =~ /\+/ ? '+' : '-');

			if ($frame =~ /\+/) {
				$frame =~ s/\+//; $frame -= 1;
			} else {
				$frame =~ s/\-//; $frame -= 1; $frame .= 'R'; # TODO: CHECK
			}

			my ($extended_lo, $extended_hi) = (0, 0);
			$extended_lo = 1 if $start < 1 or $end < 1;
			my $seq_len = length($$seqs{$cur_seq});
			$extended_hi = 1 if $start > $seq_len or $end > $seq_len;

			$start += 3 if $start < 1;
			$end += 3 if $end < 1;
			$start -= 3 if $start > $seq_len;
			$end -= 3 if $end > $seq_len;

			my $orf_len = abs($end - $start) + 1;
			die("Internal error: $start..$end [$orf_len]") if $orf_len %3 != 0;
			die("Internal error: ORF length not divisible by 3") if $orf_len %3 != 0;
#			$start = 1 if $start < 1; $end = 1 if $end < 1;
#			$start = $seq_len if $start > $seq_len; $end = $seq_len if $end > $seq_len;

			push(@{$predictions{$cur_seq}}, {seqname => $cur_seq, lo => min($start, $end), hi => max($start, $end),
			    strand => $strand, # frame => $frame,
				gl3_score => $raw_score, type => 'CDS', start => $start, stop => $end, predictor => 'Glimmer3',
				extended_lo => $extended_lo, extended_hi => $extended_hi});
		}
	}
	close PRED;
	my $total_predictions; $total_predictions += scalar(@{$predictions{$_}}) for keys %predictions;
	logmsg "Loaded $total_predictions Glimmer3 predictions in ".keys(%predictions)." sequences from file $predfile";
	return \%predictions;
}

# Apparently the default min longorf value in glimmer is 300.
# I make a wild guesstimate here and consider at least 1MB of sequence in contigs of length over 500 to be sufficient for training by default.
# This can be tuned with the parameters below.
sub checkGenePredTrainSet($$) {
	my ($seqs, $settings) = @_;
	die("Internal error: no input supplied") unless $seqs and $settings;
	$$settings{min_long_seq_length} ||= 500;
	$$settings{min_nt_in_long_seqs} ||= 1e6;

	my $total_long_seq_length;
	foreach my $seq_name (keys %$seqs) {
		next if length($$seqs{$seq_name}) < $$settings{min_long_seq_length};
		$total_long_seq_length += length($$seqs{$seq_name});
	}

	return 1 if $total_long_seq_length > $$settings{min_nt_in_long_seqs};
	return 0;
}

sub trainGenemarkModel($$) {
	my ($input_file_full, $settings) = @_;
	$$settings{tempdir} ||= AKUtils::mktempdir($settings);
	$$settings{gm_trainer} ||= AKUtils::fullPathToExec('gmsn.pl') or die;

	logmsg "Running $$settings{gm_trainer} on $input_file_full...";
	my $gm_trainer_opts = "--combine --gm $$settings{gm_trainer_xopts}";
	my $invoke_string = "cd $$settings{tempdir}; $$settings{gm_trainer} $gm_trainer_opts '$input_file_full' >/dev/null 2>&1";
	system($invoke_string);
	die("Error running $$settings{gm_trainer}: $!") if $?;
	$$settings{gmhmm_model} = "$$settings{tempdir}/GeneMark_hmm_combined.mod";
	die("Error running $$settings{gm_trainer}: no model generated") unless -f $$settings{gmhmm_model};
	return $$settings{gmhmm_model};
}

=head1 getGenemarkPredictions
Accepts filename with mfa data
Outputs GMHMMP predictions (see loadGMHMMPredictions for format)
Settings:
	gm_predictor - name of or full path to gmhmmp
	gm_trainer - name of or full path to gmsn.pl
	gm_predictor_xopts - extra options to pass to gmhmmp
	gm_trainer_xopts - extra options to pass to gmsn.pl
	gmhmm_model -  (If not specified, trains the model using GeneMarkS)
=cut
sub getGenemarkPredictions($;$) {
	my ($seqs, $settings) = @_;
	die("No input sequences supplied") if ref($seqs) ne 'HASH';
	die("gmhmm_model and gene_pred_train_file should not both be set") if $$settings{gmhmm_model} and $$settings{gene_pred_train_file};

	$$settings{tempdir} ||= AKUtils::mktempdir($settings);
	$$settings{gm_trainer} ||= AKUtils::fullPathToExec('gmsn.pl') or die;
	$$settings{gm_predictor} ||= AKUtils::fullPathToExec('gmhmmp') or die;
	$$settings{gms_datadir} ||= "/usr/share/genemarks";
	$$settings{gms_alphabet} ||= 11;
	$$settings{gm_allow_generic_heu_model} = 1 unless defined $$settings{gm_allow_generic_heu_model};

	# Save and restore sequence names
	my $i; my (%renamed_seqs, %seqname_orig2working, %seqname_working2orig);
	foreach my $seqname (keys %$seqs) {
		$seqname_orig2working{$seqname} = $seqname;
		$seqname_working2orig{$seqname} = $seqname;
		$renamed_seqs{$seqname_orig2working{$seqname}} = $$seqs{$seqname};
    next;
		$i++;
		$seqname_orig2working{$seqname} = 'seq'.$i;
		$seqname_working2orig{'seq'.$i} = $seqname;
		$renamed_seqs{$seqname_orig2working{$seqname}} = $$seqs{$seqname};
	}
	my $gm_input_file = "$$settings{tempdir}/gmhmmp_in.fasta";
	printSeqsToFile(\%renamed_seqs, $gm_input_file);

	if (not defined $$settings{gmhmm_model}) {
		if (defined $$settings{gene_pred_train_file}) {
			$$settings{gene_pred_train_file} = File::Spec->rel2abs($$settings{gene_pred_train_file});
			die("Unable to locate training data file $$settings{gene_pred_train_file}") unless -f $$settings{gene_pred_train_file};
			my $train_seqs = readMfa($$settings{gene_pred_train_file});
			die("Training data file $$settings{gene_pred_train_file} is unsuitable, probably not enough long sequences")
				unless checkGenePredTrainSet($train_seqs, $settings);

			$$settings{gmhmm_model} = trainGenemarkModel($$settings{gene_pred_train_file}, $settings);
		} elsif (checkGenePredTrainSet($seqs, $settings)) { # query data appears usable
			$$settings{gmhmm_model} = trainGenemarkModel($gm_input_file, $settings);
		} elsif ($$settings{gm_allow_generic_heu_model}) { # select precomputed model
			# also do this if model generation failed
			die("gmhmmp data directory not found in $$settings{gms_datadir}") unless -d $$settings{gms_datadir};
			my $gc_rounded = sprintf('%.0f', AKUtils::computeGC($seqs) * 100);
			$gc_rounded = 30 if $gc_rounded < 30;
			$gc_rounded = 70 if $gc_rounded > 70;
			copy($$settings{gms_datadir}.'/heuristic_mat/heu_'.$$settings{gms_alphabet}.'_'.$gc_rounded.'.mat', $$settings{tempdir})
				or die("Unable to perform copy");
			copy($$settings{gms_datadir}.'/heuristic_mod/heu_'.$$settings{gms_alphabet}.'_'.$gc_rounded.'.mod', $$settings{tempdir})
				or die("Unable to perform copy");

			$$settings{gmhmm_model} = $$settings{tempdir}.'/heu_'.$$settings{gms_alphabet}.'_'.$gc_rounded.'.mod';
			logmsg "Warning: insufficient data to train native model. Using heuristic model $$settings{gmhmm_model}";
		} else {
			die("Query data in $gm_input_file insufficient for self-training"
				. "and no generic heuristic model fallback is allowed, unable to select gmhmm model");
		}
	}
	$$settings{gmhmm_model} = File::Spec->rel2abs($$settings{gmhmm_model});
	die("Model file $$settings{gmhmm_model} not found") unless -f $$settings{gmhmm_model};

	logmsg "Running $$settings{gm_predictor} on $gm_input_file using model $$settings{gmhmm_model}...";

	my $lst_file = "$$settings{tempdir}/gm_out.lst";
  if(-f $lst_file && -s $lst_file>1){
    return loadGMHMMPredictions($lst_file, $seqs, \%seqname_working2orig);
  }
	unlink $lst_file;
	while (my ($seqname, $seq) = each(%renamed_seqs)) {
		my $temp_fna = "$$settings{tempdir}/temp.fna";
		open(OUT, '>', $temp_fna) or die("Unable to open file $temp_fna for writing: $!");
		print OUT ">$seqname\n$seq";
		close OUT;

		my $gm_predictor_opts = "-r -k -m $$settings{gmhmm_model} -o $temp_fna.lst $$settings{gm_predictor_xopts}";
    my $gmCommand="cd $$settings{tempdir}; $$settings{gm_predictor} $gm_predictor_opts $temp_fna 2>&1";
		system($gmCommand);
		die("Error running $$settings{gm_predictor} with command $gmCommand: $!") if $?;
		
		open(IN, '<', "$temp_fna.lst") or die("Internal error");
		open(LST, '>>', $lst_file);
		print LST ">$seqname\n";
		print LST while <IN>;
		close LST;
		close IN;

		unlink $temp_fna;
		unlink "$temp_fna.lst";
	}
	return loadGMHMMPredictions($lst_file, $seqs, \%seqname_working2orig);
}

=head1
Load a gmhmmp output file. Return a hash:
TODO
TODO: frame info recovery
=cut
sub loadGMHMMPredictions($;$$) {
	my ($predfile, $seqs, $seqname_remap) = @_;
	open(PRED, '<', $predfile) or die("Could not open file $predfile for reading: ".$!);
	my %predictions;
	my $cur_seq = 'unknown';
	while (<PRED>) {
		if (/^\>\s*(.+?)\s*$/) {
			$cur_seq = $1;
			$cur_seq = $$seqname_remap{$cur_seq} if(defined $seqname_remap && defined($$seqname_remap{$cur_seq}));
		} else {
			my @line = split /\s+/;
			shift(@line) if $line[0] eq '';
			next unless @line > 5 and $line[0] =~ /^\d+$/;
			my ($Gene, $Strand, $LeftEnd, $RightEnd, $GeneLength, $Class, $Spacer, $RBS_score) = @line;
			my ($extended_lo, $extended_hi) = (0, 0);
			my $frame;
			if ($LeftEnd =~ /</) {
				$extended_lo = 1;
			} elsif ($RightEnd =~ />/) {
				$extended_hi = 1;
			}
			$extended_lo = 1 if $LeftEnd =~ /</;
			$extended_hi = 1 if $RightEnd =~ />/;
			$LeftEnd =~ s/<//g; $RightEnd =~ s/>//g;

			my $orf_len = $RightEnd - $LeftEnd + 1;
			die("Internal error: ORF length not divisible by 3") if $orf_len %3 != 0;
			my ($start, $end) = ($Strand eq '+' ? ($LeftEnd, $RightEnd) : ($RightEnd, $LeftEnd));
			push(@{$predictions{$cur_seq}}, {seqname => $cur_seq, lo => $LeftEnd, hi => $RightEnd, strand => $Strand,
				pg_class => $Class, rbs_spacer => $Spacer, rbs_score => $RBS_score,
				extended_lo => $extended_lo, extended_hi => $extended_hi, 
				type => 'CDS', start => $start, stop => $end, predictor => 'gmhmmp'});
		}
	}
	close PRED;
	my $total_predictions; $total_predictions += scalar(@{$predictions{$_}}) for keys %predictions;
	logmsg "Loaded $total_predictions gmhmmp predictions in ".keys(%predictions)." sequences from file $predfile";
	return \%predictions;
}

# DEPRECATED
sub getGenePredictions($;$$$$) {
	my ($input_file, $gm_trainer, $gm_predictor, $seqs, $model) = @_;
	return getGenemarkPredictions($seqs, {gm_trainer => $gm_trainer, gm_predictor => $gm_predictor, gmhmm_model => $model});
}

sub printSeqsToFile($$;$) {
	my ($seqs, $outfile, $settings) = @_;
	open(OUT, '>', $outfile) or die "Unable to open file $outfile for writing: $!";
	my @seqnames = keys %$seqs;
	@seqnames = sort alnum @seqnames
		if $$settings{order_seqs_by_name};
	@seqnames = sort {length($$seqs{$a}) <=> length($$seqs{$b})} @seqnames
		if $$settings{order_seqs_by_length};
	foreach my $seqname (@seqnames) {
		my $seq = $$seqs{$seqname};
		$seq =~ s/(.{80})/$1\n/g;
		$seq .= "\n" unless $seq =~ /\n$/;
		print OUT ">$seqname\n$seq";
	}
	close OUT;
	return $outfile;
}

sub getNumCPUs() {
	my $num_cpus;
	open(IN, '<', '/proc/cpuinfo'); while (<IN>) { /processor\s*\:\s*\d+/ or next; $num_cpus++; } close IN;
	return $num_cpus || 1;
}

sub accessions2gb($$;$) {
	my ($accessions, $gb_filename, $settings) = @_;

	$$settings{tempdir} ||= AKUtils::mktempdir($settings);

	open(LST, '>', "$$settings{tempdir}/acc_lst") or die "Unable to open $$settings{tempdir}/acc_lst for writing: $!";
	for (my $i=0; $i<@$accessions; $i+=128) {
		my $qs = "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nucleotide&sendto=t&list_uids=";
		for (my $j=$i; $j<min($i+128, scalar(@$accessions)); $j++) {
			$qs .= $$accessions[$j].",";
		}
		print LST "$qs\n";
	}
	close LST;
	system("wget -i $$settings{tempdir}/acc_lst -O $gb_filename");
	die("Error running wget: $!") if $?;
	return $gb_filename;
}

=head1
Create a BLAST database. Return the database location.
Settings:
	formatdb_protein_db: set to true if protein sequence, false/undef if nucleotide.
	formatdb_in_place: create database in the given location, not in temporary space
=cut
sub formatBLASTdb($;$) {
	my ($input_fasta_file, $settings) = @_;
	my $input_file_full = File::Spec->rel2abs($input_fasta_file);
	my ($input_file, $input_dir) = fileparse($input_file_full);

	my $formatdb_opts = ($$settings{formatdb_protein_db} ? "-p T" : "-p F");
	$formatdb_opts .= " -t $$settings{formatdb_db_title}" if defined $$settings{formatdb_db_title};
	$formatdb_opts .= " -n $$settings{formatdb_db_name}" if defined $$settings{formatdb_db_name};
	my ($invoke_string, $db_loc);
	if ($$settings{formatdb_in_place}) {
		$invoke_string = "legacy_blast.pl formatdb $formatdb_opts -i $input_file_full";
		$db_loc = $input_file_full;
	} else {
		$$settings{tempdir} ||= AKUtils::mktempdir($settings);
		die("Directory $$settings{tempdir} does not exist") unless -d $$settings{tempdir};
		symlink($input_file_full, "$$settings{tempdir}/$input_file");
		$invoke_string = "cd $$settings{tempdir}; legacy_blast.pl formatdb $formatdb_opts -i $input_file";
		$db_loc = "$$settings{tempdir}/$input_file";
	}
	system($invoke_string); die("Error running formatdb: $!") if $?;
	return $db_loc;
}

# BLAST all sequences in %$seqs against a supplied database. Return results in tabulated format.
# mode = blastp|blastn|blastx|psitblastn|tblastn|tblastx (required)
# settings required: blast_db
# settings optional: blastseqs_xopts, num_cpus, tempdir
sub blastSeqs($$$) {
	my ($seqs, $mode, $settings) = @_;

	die("Internal error: invalid mode") if $mode !~ /^(blastp|blastn|blastx|psitblastn|tblastn|tblastx)$/;
	die("No BLAST database name supplied") unless $$settings{blast_db};

	$$settings{tempdir} ||= AKUtils::mktempdir($settings);

	$$settings{num_cpus} ||= getNumCPUs();

	my $blast_infile = "$$settings{tempdir}/blastseqs.in";
	my $blast_outfile = "$$settings{tempdir}/blastseqs.out";
	printSeqsToFile($seqs, $blast_infile);

	my $blast_qs = "legacy_blast.pl blastall -p $mode -d $$settings{blast_db} -m 8 -a $$settings{num_cpus} -i $blast_infile ";
  $blast_qs .= "-o $blast_outfile ";
	$blast_qs .= $$settings{blast_xopts};

	logmsg "Running $blast_qs" unless $$settings{quiet};
	#system($blast_qs);
  my $blsOut=`$blast_qs`;
	if ($?) {
		my $er_str = "Error running \"$blast_qs\": $!";
		$$settings{ignore_blast_errors} ? warn($er_str) : die($er_str);
	}
#	open(BLAST_OUT, "$blast_qs |") or die "Unable to run \"$blast_qs\": $!";
	return loadBLAST8($blast_outfile);
}

sub loadBLAST8($) {
	my ($blast_outfile) = @_;
	my @hits;
	open(BLAST_OUT, '<', $blast_outfile) or die "Unable to open $blast_outfile for reading: $!";
	while (<BLAST_OUT>) {
		chomp;
		my @line = split /\t/;
		my %hit;
		for (qw(name1 name2 percent_id al_len mismatch_bp gap_openings start1 end1 start2 end2 Evalue bitscore)) {
			$hit{$_} = shift @line;
		}
		push(@hits, \%hit);
	}
	close BLAST_OUT;

	return \@hits;
}

sub getBLASTGenePredictions($$) {
	my ($input_seqs, $settings) = @_;
  $$settings{num_cpus}||=$$settings{numcpus}||1;
	my $orfs = AKUtils::findOrfs2($input_seqs, {orftype=>'start2stop', need_seq=>1, need_aa_seq=>1});

  my $blastIn="$$settings{tempdir}/blastin.faa";
  my $blast_outfile="$$settings{tempdir}/blastout.bls";
  open(BLSIN,">$blastIn") or die "Could not write to $blastIn:$!";
  my $numOrfs;
  foreach my $seqname(keys %$orfs){
    foreach my $frame(keys %{$$orfs{$seqname}}){
      my $orfsInFrame=$$orfs{$seqname}{$frame};
      foreach my $stop (keys %$orfsInFrame){
        next if(!defined($$orfsInFrame{$stop})); # stops a crash if there are no genes on this contig
        my %seq=%{$$orfsInFrame{$stop}};
        # Due to the format, if there is any newlines in $seq{seq}, this format is screwed
        # Also screwed if whitespace or _ is in $seqname
        $seqname=~s/_/---/g;
        my $orfId=join("_",$seqname,$seq{start},$seq{stop},$seq{strand},$seq{ustop},$seq{length},$seq{lo},$seq{hi},$seq{seq},$seq{aa_seq});
        $orfId=~s/\s+//g;
        print BLSIN ">$orfId\n$seq{aa_seq}\n";
        $numOrfs++;
      }
    }
  }
  close BLSIN;
  logmsg "$numOrfs ORFs found.";

  # TODO limit number of results with -b -v ?
  # run with -m 9 to delimit each result with # lines
  logmsg "Running blast";
	$$settings{min_default_db_coverage} ||= 0.7;
	$$settings{min_reference_coverage} ||= 0.85;
	$$settings{ignore_blast_errors} ||= 1;
	$$settings{min_default_db_evalue} ||= 1e-8;
	$$settings{min_default_db_identity} ||= 0;
	$$settings{blast_db} ||= $$settings{prediction_blast_db};
	$$settings{blast_xopts} = " -e $$settings{min_default_db_evalue}";
	$$settings{quiet} = 1;
  # use -m 9 such that comment lines become separators
  if(!-e $blast_outfile){
    system("legacy_blast.pl blastall -p blastp -d '$$settings{blast_db}' -m 9 -a $$settings{num_cpus} -i $blastIn -o '$blast_outfile' $$settings{blast_xopts} 2>&1");
    die if $?;
    # add a comment line at the end of the blast file to mark the end of the last entry
    system("echo '# last line' >> $blast_outfile"); die if $?; 
  } else { logmsg "Found blast file; skipping blast.";}

  return parseBlastGenePredictions($blast_outfile,$settings);
}

sub parseBlastGenePredictions{
  my($blast_outfile,$settings)=@_;
  $$settings{num_cpus}||=$$settings{numcpus}||1;
  logmsg "Reading blast results";
  my $reportEvery=100;
  $reportEvery=10 if($reportEvery<10); 
  $$settings{blast_reportEvery}=$reportEvery;
  my $orfCount=0;

	my ($total, $good_cov, $rc_best); my %cov_hist;

	my %blast_preds;
  my $blastQueue=Thread::Queue->new;
  my @thread;
  push(@thread, threads->new(\&blastGenePredictionWorker,$blastQueue,$settings)) for(1..$$settings{num_cpus});
  open(BLSOUT,"<",$blast_outfile) or die "Could not open blast outfile $blast_outfile: $!";
  my @blsResult=();
  while(<BLSOUT>){
    # With -m 9, each entry is delmited by a set of comment lines.
    # Therefore queue a result whenever a comment line is found.
    if(/^#/){
      # TODO to make this safe with perl 5.8 (I think), send a string instead of an array.
      $blastQueue->enqueue(\@blsResult) if(@blsResult);
      @blsResult=();
    }else{
      push(@blsResult,$_);
    }
  }
  close BLSOUT;
  
  $blastQueue->enqueue(undef) for(@thread);
  while($blastQueue->pending>$$settings{num_cpus} && !$thread[0]->is_joinable){
    logmsg "Waiting on ".$blastQueue->pending." orfs.";
    sleep 60;
  }
  my @blast_hits=();
  for(@thread){
    my $tmp=$_->join;
    #warn "Warning: no blast hits were found in thread ".$_->tid if(!@$tmp);
    push(@blast_hits,@$tmp);
  }
  warn "WARNING: no blast hits were returned at all (internal error?)" if(!@blast_hits);
  logmsg "Found ".scalar(@blast_hits)." BLAST results";

  # @blast_hits have the best hit for every orf.
  # $blast_hits[$i]{name1} contains all the metadata.
  foreach my $best_hit (@blast_hits) {
    my($seqname,$start,$stop,$strand,$ustop,$length,$lo,$hi,$nt_seq,$aa_seq)=split(/_/,$$best_hit{name1});
    $seqname=~s/---/_/g; # change the divider back to normal
    my ($lo_in_seq, $hi_in_seq);
    if ($strand eq '+') {
      $lo_in_seq = $start + ($$best_hit{start1}-1)*3;
      $hi_in_seq = $start + ($$best_hit{end1}-1)*3;
    } else {
      $lo_in_seq = $start - ($$best_hit{end1}-1)*3;
      $hi_in_seq = $start - ($$best_hit{start1}-1)*3;
    }
    # TODO: scan up to upstream Met from hit location, but not if truncated
    my $pred = { seqname => $seqname,
      lo => ($strand eq '+' ? $lo_in_seq : $stop), # snap to stop
      hi => ($strand eq '+' ? $stop : $hi_in_seq), # snap to stop
      strand => $strand,
      blast_score => $$best_hit{bitscore},
      type => 'CDS',
      start => ($strand eq '+' ? $lo_in_seq : $hi_in_seq),
      stop => $stop,
      predictor => 'BLAST', };

    die("bad orf $$pred{lo}..$$pred{hi}") if ($$pred{hi}-$$pred{lo}+1) % 3 != 0;
    die("bad pred $$pred{lo} .. $$pred{hi} (pred $$pred{lo}..$$pred{hi}, in-hit coords $$best_hit{start1}..$$best_hit{end1})") if ($$pred{hi} - $$pred{lo} + 1 ) % 3 != 0;

    push(@{$blast_preds{$seqname}}, $pred);
  }

  # Make the blast results unique.
  # TODO figure out why they are redundant and simply prevent it.
  while(my ($seqid,$p)=each(%blast_preds)){
    my %tmpPred=();
    for(@$p){
      my $key=join("_",$$_{start},$$_{stop},$$_{strand},$$_{seqname});
      $tmpPred{$key}=$_;
    }
    $blast_preds{$seqid}=[values(%tmpPred)];
    $orfCount+=values(%tmpPred);
  }

  logmsg "Finished analyzing. Found $orfCount orfs with a good BLAST result";

	return \%blast_preds;
}

sub blastGenePredictionWorker{
  my($Q,$settings)=@_;
  #my($orfs,$Q,$settings)=@_;
  my @best_hits=();
  while(defined(my $blsResultArr=$Q->dequeue)){
    next if(!@$blsResultArr);
    chomp(@$blsResultArr);
    my @query_hits=();
    for my $result(@$blsResultArr){
      my %hit;
      my @result=split(/\s+/,$result);
      for (qw(name1 name2 percent_id al_len mismatch_bp gap_openings start1 end1 start2 end2 Evalue bitscore)) {
        $hit{$_} = shift @result;

      }
      my($seqname,$start,$stop,$strand,$ustop,$length,$lo,$hi,$nt_seq,$aa_seq)=split(/_/,$hit{name1});
      # query coverage (not db coverage)
      $hit{coverage} = $hit{al_len} * 3 / $length; 
      next if $hit{coverage} < $$settings{min_default_db_coverage};
      next if $hit{percent_id} < $$settings{min_default_db_identity} * 100;
      push(@query_hits,\%hit);
    }
    next if(!@query_hits);
    # find the best hit to add onto @best_hits. Bitscore will indicate a higher coverage and more matches.
    @query_hits=sort{$$a{bitscore} <=> $$b{bitscore}} @query_hits;
    push(@best_hits,$query_hits[0]);
  }
  return \@best_hits;
}

# Run commands using the PBS scheduler.
# Input:
#  array of jobs (required):
#   [cmd1, cmd2, ...]
#   or [{command=>cmd, [stagein=>(PBS::Client stagein format)], [stageout=>...]}, {command=>cmd, ...}, ...];
#  callback function (optional, called after each job is done with 3 arguments: job name, stdout log file, stderr log file),
#  settings hash: pbs_queue (required, name of PBS queue to use),
#   tempdir,
#   wait_on_pbs_jobs (if true, pause execution until all jobs are complete; required for callback use)
# Output: [[job1name, job1_stdout_log, job1_stderr_log], [job2name, job2_stdout_log, job2_stderr_log], ...]
# Notes:
#  Stagein is either PBS::Client stagein format or a single file or directory name.
#  Directories can be named as stagein/outs.
#  FIXME: currently requires write access to working directory, because job files are written there.
# Examples:
#  runPBSjobs(["command1", "command2"], undef, {pbs_queue => "topaz_main"});
#  runPBSjobs([{command=>"command1", stagein=>{"head:/tmp/mydir"=>"/tmp/mydir"}}],
#   $callback_subref, {pbs_queue=>"topaz_main", wait_on_pbs_jobs=>1});
#  runPBSjobs([{command=>"command1", stagein=>"/tmp/mydir", stageout=>"/tmp/mydir"}}],
#   $callback_subref, {pbs_queue=>"topaz_main", wait_on_pbs_jobs=>1});
sub runPBSjobs($$$) {
	my ($jobs, $callback, $settings) = @_;
	require PBS::Client;
	die("No commands supplied or not an array") if ref($jobs) ne 'ARRAY';
	die("No PBS queue name specified") unless defined $$settings{pbs_queue};
	die("Callback defined but is not a subroutine") if defined $callback and ref($callback) ne 'CODE';
	$$settings{tempdir} ||= AKUtils::mktempdir($settings);
	$$settings{keep_free_nodes} ||= 4;

	my %active_jobs; my @job_log;
	foreach my $i (0..$#$jobs) {
		if (ref($$jobs[$i]) ne 'HASH') {
			$$jobs[$i] = {command => $$jobs[$i]};
		}
		$$jobs[$i]->{name} ||= "runpbs.$$.$i";
	}

	logmsg "Submitting ".@$jobs." PBS jobs...";
	my $pbs_client = PBS::Client->new();
	foreach my $job (@$jobs) {
		die("Bad job name \"$$job{name}\"") if $$job{name} !~ /^[\d\w\.]+$/;
		
		if ($$settings{keep_free_nodes}) {
			while (1) {
				my $free_nodes = `pbsnodes -l free|wc -l` + 0;
				# TODO: differentiate between nodes/cpus
				# TODO: sometimes this doesnt update fast enough
				# - maybe take min($free_nodes, $initial_free_nodes - $jobs_submitted)
				last if $free_nodes > $$settings{keep_free_nodes};
				sleep 5;
			}
		}

		my $job_h = PBS::Client::Job->new(name => $$job{name},
			script => "runpbsjobs.sh",
			queue => $$settings{pbs_queue}, cmd => $$job{command},
			ppn => 8, # this makes no sense but yes, it makes it schedule one job per node at a time
			# nodes => 1,
			ofile => "$$settings{tempdir}/$$job{name}.out",
			efile => "$$settings{tempdir}/$$job{name}.err",
		);

		if (ref($$job{stagein}) =~ /(ARRAY|HASH)/) {
			$job_h->stagein($$job{stagein});
		} elsif (-e $$job{stagein}) { # FIXME: this doesnt work for directories
			$job_h->stagein({hostname().':'.$$job{stagein} => $$job{stagein}});
		} elsif (defined $$job{stagein}) {
			die("stagein \"$$job{stagein}\" defined in job $$job{name} but is not a hash, array, or filename on this machine");
		}
		if (ref($$job{stageout}) =~ /(ARRAY|HASH)/) {
			$job_h->stage($$job{stageout});
		} elsif (-e $$job{stageout}) { # FIXME: this doesnt work for directories
			$job_h->stageout({hostname().':'.$$job{stageout} => $$job{stageout}});
		} elsif (defined $$job{stageout}) {
			die("stageout \"$$job{stageout}\" defined in job $$job{name} but is not a hash, array, or filename on this machine");
		}

		my $job_id = $pbs_client->qsub($job_h);
		die("Internal error: unexpected job id @$job_id") if @$job_id != 1;
		$active_jobs{$$job{name}} = $$job_id[0].".".hostname();
		logmsg "Submitted job $$job{name} (id $$job_id[0].".hostname().")";
	}

	if ($$settings{wait_on_pbs_jobs}) {
		require PBS; # WARNING: There are two PBS.pm's out there. This uses http://www-rcf.usc.edu/~garrick/perl-PBS/
		logmsg "Waiting for PBS jobs to complete...";
		my $conn = PBS::pbs_connect(""); # TODO: specify name of pbs server
		while (%active_jobs > 0) {
			# see pbs(3) manpage for methods/symbols/ways to enum attribs
			my $status = PBS::pbs_statjob($conn, undef, [qw/job_state/], undef);
			my %busy_jobs;
			foreach my $job (@$status) {
				$busy_jobs{$$job{name}} = $$job{attribs}->[0]->{value};
			}
			foreach my $job (keys %active_jobs) {
				my $job_id = $active_jobs{$job};
				next if defined $busy_jobs{$job_id};
				$$settings{pbs_job_stageout_timeout} = 300; # in seconds
				my $t = 0;
				while (not -f "$$settings{tempdir}/$job.out") {
					$t++; die("Unable to recover stdout output for job $job") if $t > $$settings{pbs_job_stageout_timeout};
					warn "Waiting for $$settings{tempdir}/$job.out to appear\n" if $t%10 == 1;
					sleep 1;
				}
				while (not -f "$$settings{tempdir}/$job.err") {
					$t++; die("Unable to recover stderr output for job $job") if $t > $$settings{pbs_job_stageout_timeout};
					warn "Waiting for $$settings{tempdir}/$job.err to appear\n" if $t%10 == 1;
					sleep 1;
				}
				logmsg "Job $job is done"; # job is done; outprocess/callback
				if ($callback) {
					$callback->($job, "$$settings{tempdir}/$job.out", "$$settings{tempdir}/$job.err");
				}
				delete $active_jobs{$job};
				push(@job_log, [$job, "$$settings{tempdir}/$job.out", "$$settings{tempdir}/$job.err"]);
			}
			sleep 2; # TODO: subscribe to job events
		}
		logmsg "All PBS jobs done";
		PBS::pbs_disconnect($conn);
	}
	for (glob "runpbsjobs.sh*") {
		unlink; # TODO: this should probably be in END{}/sig handler block as well
	}
	return \@job_log; # FIXME: this is only populated with wait_on_pbs_jobs and contains no exit status
}

sub loadCDSFromGenbankFile($) {
    my ($file) = @_;
    my %genes;
	logmsg "Loading data from $file...";
	my $gbh = Bio::SeqIO->new(-format => 'genbank', -file => $file);
	while (my $gb_ent = $gbh->next_seq()) {
		warn "Warning: species tag isn't present in $file" unless $gb_ent->species;
		for my $cds ($gb_ent->get_SeqFeatures) {
			next if $cds->primary_tag ne 'CDS';
			my $locus_name = ($cds->get_tag_values('locus_tag'))[0];
			die("Internal error: duplicate locus id $locus_name") if $genes{$gb_ent->id}->{$locus_name};
			$genes{$gb_ent->id}->{$locus_name} = $cds;
		}
	}
	return \%genes;
}

# Load configuration values from config files.
# Usage: $settings = loadConfig($default_settings);
# By default, load 3 files in succession: bindir/appnamerc, /etc/appnamerc, ~/.appnamerc (if they exist)
# Override this with: $settings = loadConfig({config_files => [file1, file2]});
# 2011-11-22 LK added that a current working directory can have an appnamerc file too
sub loadConfig($) {
	my ($settings) = @_;

	$0 = fileparse($0) if $0 =~ /\//;
	$$settings{appname} ||= $0;
	$$settings{app_config_dir} ||= "$FindBin::RealBin/../conf";
	$$settings{system_config_dir} ||= "/etc";
	$$settings{user_config_dir} ||= $ENV{HOME};
	$$settings{config_files} ||= ["$FindBin::RealBin/$$settings{appname}rc",
		"$$settings{app_config_dir}/$$settings{appname}rc",
        "$$settings{system_config_dir}/$$settings{appname}rc",
		"$$settings{user_config_dir}/.$$settings{appname}rc",
    "./$$settings{appname}rc"];

	my @config_files = @{$$settings{config_files}}; # avoid ridiculousness from config_files being set in a config file

	foreach my $file (@config_files) {
		open(IN, '<', $file) or next;
		my $i;
		while (<IN>) {
			$i++;
			chomp;
			next if /^\s*\#/ or /^$/;
			/^\s*([^\=]+?)\s*=\s*([^\=]+?)\s*$/
				or die("Invalid configuration syntax on line $i of file $file. Expected variable=value syntax");
			my ($key, $value) = ($1, $2);
			if ($value =~ /^"/) { $value =~ s/^"//; $value =~ s/"$//; }
			$$settings{$key} = $value;
		}
		close IN;
	}

	return $settings;
}

# NB: this is inaccurate if the uncertain nucleotide designations count is a substantial fraction of total
sub computeGC($) {
	my ($seqs) = @_;
	die("Internal error: no input supplied") unless $seqs;
	my ($tot_len);
	my (%nuc_freqs);
	foreach my $seq_name (keys %$seqs) {
		$tot_len += length($$seqs{$seq_name});
		# fixme: inefficient for huge seqs
		my $seq_copy = $$seqs{$seq_name};
		foreach my $n ('A', 'T', 'G', 'C') { $nuc_freqs{$n} += ($seq_copy =~ s/$n//gi); }
	}
	return (($nuc_freqs{'G'}+$nuc_freqs{'C'})/$tot_len);
}

sub workerPool($;$) {
	my ($input_tasks, $settings) = @_;
	my $debug = 0;
	require IO::Select;

	# TODO: pass library paths properly (test by unsetting PERL5LIB)
	$$settings{worker_exec_name} ||= "perl -e'use AKUtils; exit(AKUtils::workerWrapper())'";
	die("Error: Expected an array of tasks or a hash of named tasks") if ref($input_tasks) ne 'ARRAY' and ref($input_tasks) ne 'HASH';

	my %tasks;
	if (ref($input_tasks) eq 'ARRAY') {
		my $i;
		foreach my $task (@$input_tasks) {
			$i++; $tasks{$i} = $task;
		}
	} else {
		$tasks{$_} = $$input_tasks{$_} for keys %$input_tasks;
	}

	$$settings{num_threads} ||= min(AKUtils::getNumCPUs(), scalar(keys %tasks));
	$$settings{define_thread_ids} = 1 unless defined $$settings{define_thread_ids};
	$$settings{require_report_line} = 1 unless defined $$settings{define_thread_ids};

	# Spawn the worker thread pool (not real multithreading - just processes with pipes)
	my $th_set = new IO::Select();
	my %worker_info;
	foreach my $thread_id (0..$$settings{num_threads}-1) {
		my $invoke_str = $$settings{worker_exec_name};
		if ($$settings{define_thread_ids}) {
			$$settings{thread_id_arg_name} ||= '-i';
			$invoke_str .= " $$settings{thread_id_arg_name} $thread_id";
		}
		logmsg "Launching thread $thread_id via \"$invoke_str\"";
		# open($th, '-|', $invoke_str); # $th is a read handle (worker output is piped to us)
		# open($th, '|-', $invoke_str); # $th is a write handle (we pipe input to the worker)
		use IPC::Open2;
		my ($worker_out, $worker_in);
		my $pid = open2($worker_out, $worker_in, $invoke_str);
		$worker_info{$worker_out} = {pid => $pid, rh => $worker_out, wh => $worker_in, thread_id => $thread_id};
		
		$th_set->add($worker_out);
	}

	my @outstanding_tasks = sort keys %tasks;
	my %reports;
    DISPATCH: while (my @ready_rhs = $th_set->can_read) {
		foreach my $worker_rh (@ready_rhs) {
			warn "worker rh $worker_rh ready..." if $debug;
			my $report = <$worker_rh>; chomp $report;
			$reports{$worker_info{$worker_rh}->{current_task}} = $report if defined $worker_info{$worker_rh}->{current_task};
			undef $worker_info{$worker_rh}->{current_task};
			
			# 1. accept output from worker
			# 2. if tasks remain, dispatch task, otherwise terminate this worker
			# 3. if no dispached tasks are outstanding, terminate all workers and return

			if (my $next_task_id = shift @outstanding_tasks) {
				my $next_task = $tasks{$next_task_id};
				warn "dispatching task $next_task" if $debug;
				print { $worker_info{$worker_rh}->{wh} } $next_task."\n";
				$worker_info{$worker_rh}->{current_task} = $next_task_id;
			} else {
				$th_set->remove($worker_rh);
			}
		}
	}

	logmsg "done with tasks, cleaning up\n";
	foreach my $h (keys %worker_info) {
		$th_set->remove($h);
		close $h;
		close $worker_info{$h}->{wh};
		waitpid($worker_info{$h}->{pid}, 0);
		logmsg "cleaned up worker with pid $worker_info{$h}->{pid}";
	}
	return \%reports;
}

sub workerWrapper() {
	require FileHandle;
	STDIN->autoflush(1);
	STDOUT->autoflush(1);

	warn "$0 started with args @ARGV\n";
	my $worker_id;
	$worker_id = $ARGV[1] if $ARGV[0] eq '-i';
	print "$0: worker ready\n";

	while (<STDIN>) {
		chomp;
		my $invoke_string = $_;
		system($invoke_string);
		print "$0: Command \"$invoke_string\" exited with status $?\n";
	}
	# TODO: install signal handlers
	return 0;
}

# Given an uninterrupted bacterial ORF and a position of interest,
# find the closest start codon to that position.
# Input: position of interest, sequence of interest,
# low coordinate of ORF (indexed starting at 1), high coordinate of ORF (at least one of the two is required),
# (e.g. in a sequence consisting of NATGNNNTGA, orf_lo = 2, orf_hi = 10),
# strand of ORF (+ or -) (required). Optional settings: { start_codons => 'codon1|codon2|codon3' }
# Output: the opening/most upstream coordinate (lowest for positive strand, highest for negative)
# of the start codon closest to the position of interest in the indicated frame/strand.
# The coordinate is indexed starting at 1.
sub snapToStart($$$$$;$) {
	my ($posn, $seq, $orf_lo, $orf_hi, $strand, $settings) = @_;
	if (defined $orf_lo and defined $orf_hi) {
		die "Internal error: orf_hi < orf_lo" if $orf_hi < $orf_lo;
		die "Internal error: invalid orf coordinates" if (($orf_hi - $orf_lo + 1) % 3 != 0);
	}
	die "Internal error: No strand argument supplied" unless $strand eq '+' or $strand eq '-';
	die "Internal error: ORF coordinate exceeds sequence length" if defined $orf_lo and length($seq) < $orf_lo;
	die "Internal error: ORF coordinate exceeds sequence length" if defined $orf_hi and length($seq) < $orf_hi;
	die "Internal error: ORF coordinate less than 1" if defined $orf_hi and $orf_hi < 1;
	die "Internal error: ORF coordinate less than 1" if defined $orf_lo and $orf_lo < 1;

	$$settings{start_codons} ||= "ATG|GTG";
	my $start_codons;
	if ($strand eq '+') {
		$start_codons = $$settings{start_codons};
	} else {
		$start_codons = reverse($$settings{start_codons}); $start_codons =~ tr/ATGC/TACG/;
	}

	my $frame; # NB: no 'R' suffix in this block
	if (defined $orf_lo) {
		$frame = ($orf_lo - 1) % 3;
		# $frame .= 'R' if $strand eq '-';
	} elsif (defined $orf_hi) {
		$frame = $orf_hi % 3;
		# $frame .= 'R' if $strand eq '-';
	}
	$orf_lo ||= $frame+1;
	$orf_hi ||= length($seq); # FIXME - length($seq) % 3;
	
	my ($prev_start_lo, $next_start_lo);
	for (my $next_codon_lo = $posn-1 - (($posn-$frame-1) % 3); $next_codon_lo < $orf_hi; $next_codon_lo += 3) {
		next if $next_codon_lo < 0;
		my $codon = substr($seq, $next_codon_lo, 3);
		last if length($codon) < 3;
		print "frame $frame @ pos=$next_codon_lo Next codon $codon\n";
		if ($codon =~ /($start_codons)/) {
			warn "found start @ $next_codon_lo\n";
			$next_start_lo = $next_codon_lo; last;
		}
	}
	for (my $prev_codon_lo = $posn-1 - (($posn-$frame-1) % 3) - 3; $prev_codon_lo > $orf_lo-1; $prev_codon_lo -= 3) {
		my $codon = substr($seq, $prev_codon_lo, 3);
		last if length($codon) < 3;
		print "frame $frame @ pos=$prev_codon_lo Prev codon $codon\n";
		if ($codon =~ /($start_codons)/) {
			warn "found start @ $prev_codon_lo\n";
			$prev_start_lo = $prev_codon_lo; last;
		}
	}
	warn "start codons $start_codons; psl $prev_start_lo  nsl $next_start_lo\n";
	my $best_start;
	if (defined $prev_start_lo and defined $next_start_lo) {
		if (($posn-1) - $prev_start_lo < $next_start_lo - ($posn-1)) {
			$best_start = $prev_start_lo;
		} else {
			$best_start = $next_start_lo;
		}
	} else {
		$best_start = $next_start_lo if defined $next_start_lo;
		$best_start = $prev_start_lo if defined $prev_start_lo;
	}

	$best_start += 1 if defined $best_start; # convert from 0 to 1 based
	$best_start += 2 if $strand eq '-' and defined $best_start;

	return $best_start;
}

#################
## LK subs
#################
# See whether a fastq file is paired end or not. It must be in a velvet-style shuffled format.
# In other words, the left and right sides of a pair follow each other in the file.
# params: fastq file and settings
# fastq file can be gzip'd
# settings:  checkFirst is an integer to check the first X deflines
# TODO just extract IDs and send them to the other _sub() 
sub is_fastqPE($;$){
  my($fastq,$settings)=@_;

  # if checkFirst is undef or 0, this will cause it to check at least the first 20 entries.
  $$settings{checkFirst}||=20;
  $$settings{checkFirst}=20 if($$settings{checkFirst}<2);

  # it is paired end if it validates with any naming system
  my $is_pairedEnd=_is_fastqPESra($fastq,$settings) || _is_fastqPECasava18($fastq,$settings) || _is_fastqPECasava17($fastq,$settings);
  
  return $is_pairedEnd;
}

sub _is_fastqPESra{
  my($fastq,$settings)=@_;
  my $numEntriesToCheck=$$settings{checkFirst}||20;
  
  my $numEntries=0;
  my $fp;
  if($fastq=~/\.gz$/){
    open($fp,"gunzip -c '$fastq' |") or die "Could not open $fastq for reading: $!";
  }else{
    open($fp,"<",$fastq) or die "Could not open $fastq for reading: $!";
  }
  my $discard;
  while(<$fp>){
    chomp;
    s/^@//;
    my($genome,$info1,$info2)=split(/\s+/);
    if(!$info2){
      close $fp;
      return 0;
    }
    my($instrument,$flowcellid,$lane,$x,$y,$X,$Y)=split(/:/,$info1);
    $discard=<$fp> for(1..3); # discard the sequence and quality of the read for these purposes
    my $secondId=<$fp>;
    my($genome2,$info3,$info4)=split(/\s+/,$secondId);
    my($instrument2,$flowcellid2,$lane2,$x2,$y2,$X2,$Y2)=split(/:/,$info3);
    $_||="" for($X,$Y,$X2,$Y2); # these variables might not be present
    if($instrument ne $instrument2 || $flowcellid ne $flowcellid2 || $lane ne $lane2 || $x ne $x2 || $y ne $y2 || $X ne $X2 || $Y ne $Y2){
      close $fp;
      return 0;
    }
    $discard=<$fp> for(1..3);
    $numEntries+=2;
    last if($numEntries > $numEntriesToCheck);
  }
  return 1;
}

sub _is_fastqPECasava18{
  my($fastq,$settings)=@_;
  my $numEntriesToCheck=$$settings{checkFirst}||20;
  
  my $numEntries=0;
  my $fp;
  if($fastq=~/\.gz$/){
    open($fp,"gunzip -c '$fastq' |") or die "Could not open $fastq for reading: $!";
  }else{
    open($fp,"<",$fastq) or die "Could not open $fastq for reading: $!";
  }
  while(<$fp>){
    chomp;
    s/^@//;
    my($instrument,$runid,$flowcellid,$lane,$tile,$x,$yandmember,$is_failedRead,$controlBits,$indexSequence)=split(/:/,$_);
    my $discard;
    $discard=<$fp> for(1..3); # discard the sequence and quality of the read for these purposes
    my($y,$member)=split(/\s+/,$yandmember);
    
    # if all information is the same, except the member (1 to 2), then it is still paired until this point.
    my $secondId=<$fp>;
    chomp $secondId;
    $secondId=~s/^@//;
    my($inst2,$runid2,$fcid2,$lane2,$tile2,$x2,$yandmember2,$is_failedRead2,$controlBits2,$indexSequence2)=split(/:/,$secondId);
    $discard=<$fp> for(1..3); # discard the sequence and quality of the read for these purposes
    my($y2,$member2)=split(/\s+/,$yandmember2);

    if($instrument ne $inst2 || $runid ne $runid2 || $flowcellid ne $fcid2 || $tile ne $tile2 || $member!=1 || $member2!=2){
      #logmsg "Failed!\n$instrument,$runid,$flowcellid,$lane,$tile,$x,$yandmember,$is_failedRead,$controlBits,$indexSequence\n$inst2,$runid2,$fcid2,$lane2,$tile2,$x2,$yandmember2,$is_failedRead2,$controlBits2,$indexSequence2\n";
      close $fp;
      return 0;
    }

    $numEntries+=2;
    last if($numEntries>$numEntriesToCheck);
  }

  close $fp;

  return 1;
}

sub _is_fastqPECasava17{
  my($fastq,$settings)=@_;
  # 20 reads is probably enough to make sure that it's shuffled (1/2^20 chance I'm wrong)
  my $numEntriesToCheck=$$settings{checkFirst}||20;
  my $numEntries=0;
  my $fp;
  if($fastq=~/\.gz$/){
    open($fp,"gunzip -c '$fastq' |") or die "Could not open $fastq for reading: $!";
  }else{
    open($fp,"<",$fastq) or die "Could not open $fastq for reading: $!";
  }
  while(my $read1Id=<$fp>){
    my $discard;
    $discard=<$fp> for(1..3);
    my $read2Id=<$fp>;
    $discard=<$fp> for(1..3);

    if($read1Id!~/\/1$/ || $read2Id!~/\/2$/){
      close $fp;
      return 0;
    }

    $numEntries+=2;
    last if($numEntries>=$numEntriesToCheck);
  }
  close $fp;

  return 1;
}


__END__

# + strand only for now
# seq coords start at 1
sub startStopPos($) {
	my ($seq) = @_;
	my ($start_codons, $stop_codons) = (['ATG','GTG'], ['TAA','TAG','TGA']);
	my (%start_codon_pos, %stop_codon_pos);

	foreach my $pos (0..length($seq)/3) {
		FRAME: foreach my $frame (0, 1, 2) {
			my $codon = substr($seq, ($pos*3)+$frame, 3);
			foreach my $startc (@$start_codons) {
				if ($codon eq $startc) {
					push(@{$start_codon_pos{$frame}}, ($pos*3)+$frame+1);
					next FRAME;
				}
			}
			foreach my $stopc (@$stop_codons) {
				if ($codon eq $stopc) {
					push(@{$stop_codon_pos{$frame}}, ($pos*3)+$frame+1);
					next FRAME;
				}
			}
		}
	}
	return (\%start_codon_pos, \%stop_codon_pos);
}
