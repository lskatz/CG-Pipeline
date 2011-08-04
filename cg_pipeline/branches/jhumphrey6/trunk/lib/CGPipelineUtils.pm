#!/usr/bin/env perl

# AssemblyUtils.pm: Shotgun genome assembly-related utility subroutines
# Author: Andrey Kislyuk (kislyuk@gatech.edu)

package AssemblyUtils;
require 5.005;
my ($VERSION) = ('$Id$' =~ /,v\s+(\d+\S+)/o);

use strict;
use Exporter;
use IPC::Open2;
use Cwd;
use List::Util qw(min max sum reduce shuffle);
use File::Basename;
use File::Spec;
use File::Temp ('tempdir');
use Sys::Hostname;
use Storable ('dclone'); # for deep copying

local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };

our @ISA = "Exporter";
our @methods = qw();
our %EXPORT_TAGS = (all => [@methods]);
Exporter::export_ok_tags('all');

sub logmsg {my $FH = $AKUtils::LOG || *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

1;

