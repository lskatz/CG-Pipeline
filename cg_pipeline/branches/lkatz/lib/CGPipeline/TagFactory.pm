#!/usr/bin/env perl

# CGPipeline::TagFactory: Generate genomic feature tags for gene prediction and annotation of genomes.
# Author: Andrey Kislyuk (kislyuk@gatech.edu)

package CGPipeline::TagFactory;
require 5.005;
my ($VERSION) = ('$Id$' =~ /,v\s+(\d+\S+)/o);

use strict;
use Exporter;
use List::Util qw(min max sum reduce shuffle);

local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };

our @ISA = "Exporter";
our @methods = qw();
our %EXPORT_TAGS = (all => [@methods]);
Exporter::export_ok_tags('all');

sub logmsg {my $FH = $AKUtils::LOG || *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

sub new {
	my ($class, $fields) = @_;
	$class = ref($class) || $class;
	my $self = $fields;
	die("field strain_name is required") unless $$self{strain_name};
	die("field factory_type is required") unless $$self{factory_type};

	die("Unknown factory type $$self{factory_type}. Implemented factory types: draft_orf_tagger")
		if $$self{factory_type} ne "draft_orf_tagger";

	$$self{tag_prefix} ||= $$self{strain_name}."_";
	$$self{idx_counter} = 0;

    bless $self, $class;
    return $self;
}

sub nextTag {
	my ($self) = @_;

	$$self{idx_counter}++;
	
	return $$self{tag_prefix}.sprintf('%.4d', $$self{idx_counter});
}


1;
