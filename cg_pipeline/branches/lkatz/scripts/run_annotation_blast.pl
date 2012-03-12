#!/usr/bin/env perl

# run-assembly: Perform standard assembly protocol operations on 454 pyrosequenced flowgram file(s)
# Author: Andrey Kislyuk (kislyuk@gatech.edu)
# Author: Lee Katz <lkatz@cdc.gov>

package PipelineRunner;
my ($VERSION) = ('$Id: $' =~ /,v\s+(\d+\S+)/o);

my $settings = {
  appname => 'cgpipeline',
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

  my @cmd_options = qw(outfile=s blastfile=s db=s parametersForBlast=s min_aa_coverage=i min_aa_identity=i min_aa_similarity=i tempdir=s keep);
  GetOptions($settings, @cmd_options) or die;

  $$settings{blast_db}||=File::Spec->rel2abs($$settings{db});

  die(usage()) if @ARGV != 1;
  $$settings{min_aa_coverage}||=1;
  $$settings{min_aa_identity}||=1;
  $$settings{min_aa_similarity}||=1;
  for (qw(blast_db)) {
    die("Argument $_ must be supplied") unless $$settings{$_};
  }
  $$settings{query_mfa} = $ARGV[0];
  die("File $$settings{query_mfa} does not exist") unless -f $$settings{query_mfa};
  $$settings{outfile} ||= "$$settings{query_mfa}.unspecifiedhits.sql";

  $$settings{tempdir} ||= tempdir(File::Spec->tmpdir()."/$0.$$.XXXXX", CLEANUP => !($$settings{keep}));
  logmsg "Temporary directory is $$settings{tempdir}";
  my $numcpus=AKUtils::getNumCPUs();

  # don't blast if a blastfile has been given
  my $report;
  if(!$$settings{blastfile}){
    $ENV{BLASTDB} = (fileparse($$settings{blast_db}))[1];
    my $command="legacy_blast.pl blastall -p blastp -a $numcpus -o $$settings{tempdir}/$0.$$.blast_out -d $$settings{blast_db} -i $$settings{query_mfa} $$settings{parametersForBlast} -m 0";
    logmsg "Running BLAST on $$settings{query_mfa} vs. $$settings{blast_db}...\n  $command";
    system($command);
    die "Problem with blast" if $?;
    logmsg "Finished with BLAST";
    $report=new Bio::SearchIO(-format=>'blast',-file=>"$$settings{tempdir}/$0.$$.blast_out");
  } else {
    $report=new Bio::SearchIO(-format=>'blast',-file=>$$settings{blastfile});
  }

  my $resultCounter=0;
  logmsg "Reading $$settings{query_mfa} and parsing blast result";
  my %query_seqs;
  #my %query_seqs :shared;
  my $query_seqs = AKUtils::readMfa($$settings{query_mfa}); # should be shared
  %query_seqs=%$query_seqs;

  my %reported_hits;

  # launch threads for reading the output file
  my $Q=Thread::Queue->new;
  my $printQueue=Thread::Queue->new;
  my @thr;
  $thr[$_]=threads->new(\&readBlastOutputWorker,$Q,$printQueue,\%query_seqs,$settings) for(0..$numcpus-1);
  my $printWorkerThread=threads->new(\&printWorker,$printQueue,$settings); # prints the sequences to file

  while (my $result = $report->next_result) {
    $resultCounter++; 
    while (my $hit = $result->next_hit) {
      while (my $hsp = $hit->next_hsp) {
        eval{
          my($i1,$i2)=($hsp->seq('hit')->id, $hsp->seq('query')->id);
        };
        if($@){
          logmsg "Warning: problem with an hsp. Moving on.\n";
          next;
        }
        my $reported_hit={
          query_id => $hsp->seq('query')->id, 
          target_id => $hsp->seq('hit')->id, 
          evalue => $hsp->evalue,
          db_name => $$settings{blast_db},
          percent_identity => $hsp->percent_identity,
          description=>$hit->description || '.',
          rank=>$hit->rank,
          score=>$hit->score,
          bits=>$hit->bits,
          percent_conserved=>$hsp->frac_conserved*100,
          hitName=>$hit->name,
          hspLength=>$hit->length,
          hspQueryString=>$hsp->query_string,
          hspHitString=>$hsp->hit_string,
          hspHomologyString=>$hsp->homology_string,
        };
        $Q->enqueue($reported_hit);
      }
    }
    logmsg "Processed $resultCounter proteins" if($resultCounter % 100 == 0);
    last if $resultCounter>20;
  }
  $Q->enqueue(undef) for(0..$numcpus-1);
  $printQueue->enqueue(undef);

  # status update on the queues
  while($Q->pending || $printQueue->pending){
    my($qPending,$pPending)=($Q->pending,$printQueue->pending);
    logmsg "Reading blast queue pending: $qPending; Printing results pending: $pPending";
    sleep 5;
    if($Q->pending == $qPending && $printQueue->pending == $pPending){
      logmsg "Warning: There hasn't been much progress with the queues. I'm finishing these queue status updates in case it hurries this script along.";
      last;
    }
  }

  logmsg "Joining threads";
  $_->join for(@thr);
  $printWorkerThread->join;

  logmsg "Report is in $$settings{outfile}";
  return 0;
}

#############
#############
sub readBlastOutputWorker{
  my($Q,$printQueue,$query_seqs,$settings)=@_;
  while(defined(my $r=$Q->dequeue)){ # $r is $reported_hit
    (undef, $$r{hit_accession}, $$r{hit_name}) = split(/\|/, $$r{hitName});

    my $queryHsp=$$r{hspQueryString};
    $queryHsp=~s/\-|\s//g;
    my $l2=length($$query_seqs{$$r{query_id}}); die("Internal error - something wrong with the query seq") unless $l2;
    $$r{query_coverage}=length($queryHsp)/$l2*100;
    $$r{query_coverage}=sprintf("%.2f", $$r{query_coverage});
    #$$r{query_coverage}=length($$r{hspQueryString})/$l2*100; 
    #$$r{query_coverage}=100 if($$r{query_coverage}>100);
    $$r{percent_identity} = sprintf("%.2f", $$r{percent_identity});

    # sanity checks
    if($$r{target_id} eq $$r{query_id}){
      warn("Internal error - $$r{target_id} hit against itself. I will not parse this result.");
      next;
    }
    if ($$r{query_coverage} > 100){
      warn("Internal error - coverage is > 100% for the HSP between $$r{target_id} and $$r{query_id} (Coverage: $$r{query_coverage})\nThe HSP:\n".join("\n",sprintf("%5.5s ",$$r{target_id}).$$r{hspHitString},"      ".$$r{hspHomologyString},sprintf("%5.5s ",$$r{query_id}).$$r{hspQueryString},""));
      next;
    }

    if ($$r{query_coverage} > $$settings{min_aa_coverage}
      and $$r{percent_conserved} >= $$settings{min_aa_similarity}
      and $$r{percent_identity}  >= $$settings{min_aa_identity}) {

      # send this reported hit to a second Queue to write the results
      $printQueue->enqueue($r);
    }
  }
  return 1;
}

sub printWorker{
  my($printQueue,$settings)=@_;

  my $escapeCharacter="escapeXYZescape";
  my $tmpOutfile="$$settings{tempdir}/unsorted.sql";
  logmsg "Writing to $tmpOutfile";
  open(OUT, '>', $tmpOutfile) or die;
  while(defined(my $hit=$printQueue->dequeue)){
    my @l;
    push(@l, $$hit{$_}) for qw(query_id target_id evalue query_coverage db_name percent_identity hspLength description rank score bits percent_conserved hit_accession hit_name);
    s/\|/$escapeCharacter/g for @l; # escape the pipe characters
    print OUT join('|', @l)."\n";
  }
  close OUT;

  logmsg "Sorting results to final output, $$settings{outfile}";
  system("sort -t '|' -k 1,1 -k 3,3n $tmpOutfile > $$settings{outfile}");
  #system("sort -k 1,3 -n $tmpOutfile | sed 's/$escapeCharacter/\\\\|/g' > $$settings{outfile}");
  die "Problem with sort/sed: $!" if $?;

  return 1;
}


sub usage{
  "Usage: $0 input.mfa
  -d blastdatabase
    blast formatted database
  -p 'blast parameters'
    A quoted string to pass to blastall (legacy_blast.pl blastall -p blastp is used internally)
    There is no sanity check for this parameter.
  -o outfile.sql
    pipe delimited output file

  -b blast results file
    Optionally to skip blast and use your own blast file. In human-readable form (-m 0)
  --min_aa_coverage integer 1 to 100
  --min_aa_identity integer 1 to 100
  --min_aa_similarity integer 1 to 100
  -t tempdir (default is a subdir in /tmp)
  "
}
