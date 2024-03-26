#!/usr/bin/env perl

# run-annotation-phast: run blast against a phage database and determine whether there is a phage present
# Author: Lee Katz <lkatz@cdc.gov>

package PipelineRunner;
my ($VERSION) = ('$Id: $' =~ /,v\s+(\d+\S+)/o);

my $settings = {
  appname => 'cgpipeline',
  escapeCharacter=>"escapeXYZescape", # for regular expressions
};
my $stats;

use strict;
use FindBin;
use lib "$FindBin::RealBin/../lib";
$ENV{PATH} = "$FindBin::RealBin:".$ENV{PATH};
use AKUtils qw(logmsg);

use Getopt::Long;
use File::Temp ('tempdir');
use File::Path;
use File::Spec;
use File::Copy;
use File::Basename;
use List::Util qw(min max sum shuffle);
use CGPipelineUtils;
use Bio::SearchIO;
use Data::Dumper;

use threads;
use threads::shared;
use Thread::Queue;

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };

exit(main());

sub main() {
  $settings = AKUtils::loadConfig($settings);

  die usage() if(@ARGV<1);
  my @cmd_options = qw(outfile=s db=s tempdir=s keep help numcpus=i);
  GetOptions($settings, @cmd_options) or die;
  die usage() if($$settings{help});
  $$settings{numcpus}||=1;
  $$settings{db}||=$$settings{phast_blast_db};
  $$settings{db}||=File::Spec->rel2abs($$settings{db});

  die("ERROR: cannot find a blast database at $$settings{db}. Set it using phast_blast_db in the config file (cgpipelinerc) or by using the -d setting for this script. If being run from run_annotation, it is possible that the database has not been set for this specific task.\n".usage()) if(!-e "$$settings{db}.pin" && !-e "$$settings{db}.pal" && !-e "$$settings{db}.pal");

  die("ERROR: ARGV!=1: ".join(" ",@ARGV)."\n".usage()) if @ARGV != 1;
  for (qw(db)) {
    die("ERROR: Argument $_ must be supplied") unless $$settings{$_};
  }
  $$settings{query_mfa} = $ARGV[0];
  die("File $$settings{query_mfa} does not exist") unless -e $$settings{query_mfa};
  $$settings{outfile} ||= "$$settings{query_mfa}.phast.sql";

  $$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($$settings{keep}));
  logmsg "Temporary directory is $$settings{tempdir}";

  my $genePreds=findPhastGenes($settings);
  # see if there are consecutive phage genes from the same phage

  return 0;
}

sub findPhastGenes{
  my($settings)=@_;
  
  my $blsFile="$$settings{tempdir}/phast.asm.bls";
  if(!-e $blsFile){
    my $command="blastx -db $$settings{db} -query $$settings{query_mfa} -num_threads $$settings{numcpus} -evalue 0.05 -max_target_seqs 1 -outfmt 6 > $blsFile";
    logmsg $command;
    system($command);
    die if $?;
  }else{
    logmsg "Found $blsFile; not recreating";
  }

  my %phageLoc=();
  my @header=qw(qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore);
  open(BLS,$blsFile) or die "ERROR: Could not open $blsFile:$!";
  while(<BLS>){
    chomp;
    my %F;
    @F{@header}=split /\t/,$_;

    # phage=>{start}={...}
    $phageLoc{$F{qseqid}}{$F{qseqid}}{$F{qstart}}=\%F;
  }
  close BLS;
  die Dumper \%phageLoc;
  return \%phageLoc;

}

sub usage{
  "Usage: $0 input.assembly.fasta
  -d blastdatabase
    protein blast formatted database
  -o outfile.sql
    pipe delimited output file
  --numcpus 1 The number of processors to use
  -t tempdir
  "
}
