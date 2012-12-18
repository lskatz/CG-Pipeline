#!/usr/bin/env perl

my $settings = {};
my $stats = {};

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
use Bio::SeqIO;
use XML::LibXML::Reader;

use threads;
use Thread::Queue;

$0 = fileparse($0);

exit(main());

sub main{
  my $usage = "$0 [-outfile=parsed_interpro_output.sql] input.mfa";
  my @cmd_options = ('tempdir=s', 'outfile=s', 'keep');
  GetOptions($settings, @cmd_options) or die "Usage: $usage";
  die("Usage: $usage") if @ARGV != 1;
  #die("Argument blast_db must be supplied") unless $$settings{blast_db};
  $$settings{query_mfa} = $ARGV[0];
  die("File $$settings{query_mfa} does not exist") unless -f $$settings{query_mfa};
  $$settings{outfile} ||= "$$settings{query_mfa}.iprscan_out.xml";

  my $updateQueue=Thread::Queue->new;
  my $updateThr=threads->create(\&monitorOutfile,$$settings{outfile},$updateQueue,$settings);

  my $invoke_string = "iprscan -cli -i '$$settings{query_mfa}' -o '$$settings{outfile}' -iprlookup -goterms";
  $invoke_string .= " -altjobs";
  logmsg "Running \"$invoke_string\"...";
  system($invoke_string);
  die if $?;

  logmsg "Parsing iprscan results in $$settings{outfile}...\n";
  my $invoke_string = "run_annotation_interpro_parse_xml.pl '$$settings{outfile}'";
  system($invoke_string);
  die if $?;

  # send a term signal to the queue
  $updateQueue->enqueue(undef);
  $updateThr->join;
  return 0;
}

sub monitorOutfile{
  my($outfile,$Q,$settings)=@_;
  sleep 5 while(!-e $outfile);
  while($Q->pending<1){
    my $numFinished=`grep -c "<protein" '$outfile'`+0;
    logmsg "Scanned approximately $numFinished proteins so far";
    for(1..30){
      last if($Q->pending>0);
      sleep 2;
    }
  }
  logmsg "Finished with interproscan";
  return 1;
}
