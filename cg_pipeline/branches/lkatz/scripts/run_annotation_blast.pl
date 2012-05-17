#!/usr/bin/env perl

# run-assembly: Perform standard assembly protocol operations on 454 pyrosequenced flowgram file(s)
# Author: Andrey Kislyuk (kislyuk@gatech.edu)
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
use Bio::SearchIO::Writer::TextResultWriter;
use Bio::SearchIO::Writer::HTMLResultWriter;

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

  $$settings{blast_db}=$$settings{db} if($$settings{db});
  $$settings{blast_db}||=File::Spec->rel2abs($$settings{db});
  die("ERROR: cannot find a blast database at $$settings{blast_db}. Set it using blast_db in the config file (cgpipelinerc) or by using the -d setting for this script. If being run from run_annotation, it is possible that the database has not been set for this specific task.\n".usage()) if(!-e "$$settings{blast_db}.pin" && !-e "$$settings{blast_db}.pal" && !-e "$$settings{blast_db}.pal");

  die("ERROR: ARGV!=1: ".join(" ",@ARGV)."\n".usage()) if @ARGV != 1;
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
  #$numcpus=1; logmsg "DEBUGGING NUMBER OF CPUS";

  # don't blast if a blastfile has been given
  my $report;
  if(!$$settings{blastfile}){
    # Progress report for the blast that will happen
    my $progressQ=Thread::Queue->new;
    my $progressThread=threads->new(\&blastProgressUpdater,$progressQ,$settings);

    $ENV{BLASTDB} = (fileparse($$settings{blast_db}))[1];
    my $command="legacy_blast.pl blastall -p blastp -a $numcpus -o $$settings{tempdir}/$0.$$.blast_out -d $$settings{blast_db} -i $$settings{query_mfa} $$settings{parametersForBlast} -m 0";
    logmsg "Running BLAST on $$settings{query_mfa} vs. $$settings{blast_db}...\n  $command";
    system($command);
    die "Problem with blast" if $?;
    $progressQ->enqueue(undef); # send term signal to progress thread
    $progressThread->join;

    $report=new Bio::SearchIO(-format=>'blast',-file=>"$$settings{tempdir}/$0.$$.blast_out");
    $$settings{blastfile}="$$settings{tempdir}/$0.$$.blast_out";
  } else {
    $report=new Bio::SearchIO(-format=>'blast',-file=>$$settings{blastfile});
  }

  logmsg "Reading $$settings{query_mfa} and parsing blast result, using $numcpus threads";
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

  my $resultCounter=0;
  local $/="\nQuery=";
  open(BLASTOUT,$$settings{blastfile}) or die "Could not open the blast output file for reading: $!";
  <BLASTOUT>; # discard the first chunk because it is not a result
  while(my $entry=<BLASTOUT>){
    $resultCounter++; 
    $entry="Query=$entry"; #put back in the query=
    my $entryFile="$$settings{tempdir}/result$resultCounter.$$.bls";
    open(ENTRY,">$entryFile") or die "Could not write to file $entryFile: $!";
    print ENTRY $entry;
    close ENTRY;

    $Q->enqueue($entryFile);
    logmsg "Enqueued $resultCounter proteins. ".$Q->pending." awaiting processing. ".$printQueue->pending." in the print queue." if($resultCounter % 100 == 0);
  }
  close BLASTOUT;

  # status update on the queues
  while($Q->pending || $printQueue->pending){
    my($qPending,$pPending)=($Q->pending,$printQueue->pending);
    logmsg "Reading blast queue pending: $qPending; Printing results pending: $pPending";
    sleep 5;
  }

  # send termination signal to threads
  $Q->enqueue(undef) for(0..$numcpus-1);
  logmsg "Joining threads";
  $_->join for(@thr);

  $printQueue->enqueue(undef); # term signal to print, now that other threads are done
  logmsg "Joining print worker thread";
  $printWorkerThread->join;

  logmsg "Report is in $$settings{outfile}";
  return 0;
}

#############
#############
sub readBlastOutputWorker{
  my($Q,$printQueue,$query_seqs,$settings)=@_;
  while(defined(my $resultFile=$Q->dequeue)){ 
    # This instance of SearchIO does not read a blast output file with a header.
    # I don't THINK it's a problem but I could be wrong.
    # If it is a problem, then the header information is found in $$settings{blastHeader}.
    # However, concatenating the header with the infile is somehow time intensive and therefore isn't performed.
    my $searchIn=Bio::SearchIO->new(-file=>$resultFile,-format=>"blast");

    while(my $result=$searchIn->next_result){

      while (my $hit = $result->next_hit) {
        while (my $hsp = $hit->next_hsp) {
          #my($query_id,$target_id,$hspQueryString
          my %r=(
            query_id => $hsp->seq('query')->id,         # printed
            target_id => $hsp->seq('hit')->id,          # printed
            hspQueryString=>$hsp->query_string,         # used for coverage
          );
          # sanity check filter
          next if($r{target_id} eq $r{query_id});
          # coverage filter
          my $queryHsp=$r{hspQueryString};
          $queryHsp=~s/\-|\s//g;
          my $l2=length($$query_seqs{$r{query_id}}); die("Internal error - something wrong with the query seq") unless $l2;
          $r{query_coverage}=length($queryHsp)/$l2*100;
          $r{query_coverage}=sprintf("%.2f", $r{query_coverage});
          next if($r{query_coverage} < $$settings{min_aa_coverage});
          # percent_conserved filter
          $r{percent_conserved}=$hsp->frac_conserved*100;
          next if($r{percent_conserved} < $$settings{min_aa_similarity});
          # percent_identity filter
          $r{percent_identity}=$hsp->percent_identity;
          next if($r{percent_identity}  < $$settings{min_aa_identity});

          # PASSED: load up the rest of the variables for analysis
          %r=(
            evalue => $hsp->evalue,                     # printed
            db_name => $$settings{blast_db},            # printed
            description=>$hit->description || '.',      # printed
            rank=>$hit->rank,                           # printed
            score=>$hit->score,                         # printed
            bits=>$hit->bits,                           # printed
            hitName=>$hit->name,                        # used for hit_accession and hit_name
            hspLength=>$hit->length,                    # printed
            %r,
          );
          (undef, $r{hit_accession}, $r{hit_name}) = split(/\|/, $r{hitName});
          $r{hit_accession}||=$r{hitName};
          $r{hit_name}||=$r{hit_accession};
          $r{percent_identity} = sprintf("%.2f", $r{percent_identity});

          # sanity checks
          if ($r{query_coverage} > 100){
            warn("Internal error - coverage is > 100%. Skipping this HSP between $r{query_id} and $r{target_id}\n");
            next;
          }

          # send this reported hit to a second Queue to write the results
          my $escapeCharacter=$$settings{escapeCharacter};
          my @l;
          push(@l, $r{$_}) for qw(query_id target_id evalue query_coverage db_name percent_identity hspLength description rank score bits percent_conserved hit_accession hit_name);
          s/\|/$escapeCharacter/g for @l; # escape the pipe characters
          my $printLine=join("|", @l)."\n";
          $printQueue->enqueue($printLine);
        }
      } # END while($hit)
    }   # END while($result)
    close BLASTIN;
  }     # END file=$q->dequeue
  return 1;
}

# receives things to print; prints them.
# Then, sorts the file better using GNU/Linux
sub printWorker{
  my($printQueue,$settings)=@_;

  my $escapeCharacter=$$settings{escapeCharacter};
  my $tmpOutfile="$$settings{tempdir}/unsorted.sql";
  logmsg "Writing to temp unsorted file $tmpOutfile";
  open(OUT, '>', $tmpOutfile) or die;
  while(defined(my $printLine=$printQueue->dequeue)){
    print OUT $printLine;
  }
  close OUT;

  logmsg "Sorting results to final output, $$settings{outfile}";
  system("sort -t '|' -k 1,1 -k 3,3g $tmpOutfile | sed 's/$escapeCharacter/\\\\|/g' > $$settings{outfile}");
  die "Problem with sort/sed: $!" if $?;

  unlink($tmpOutfile);
  warn("Warning: could not remove temp file $tmpOutfile because $!") if $?;

  return 1;
}

sub blastProgressUpdater{
  my($Q,$settings)=@_;
  my $outfile="$$settings{tempdir}/$0.$$.blast_out";
  # wait until a file is there before updating
  sleep 5 while(!-e $outfile);
  sleep 5;
  while($Q->pending==0){ # the queue will have one item in it to signal for it to be done
    my $numFinished=`grep -c "Query=" '$outfile'`+0;
    logmsg "Finished with $numFinished queries";
    # give it a chance to break out of this subroutine if it has been given the term signal
    # Hold for 1 minute before updating
    for(1..30){
      last if($Q->pending>0);
      sleep 2;
    }
  }
  logmsg "Finished with BLAST";
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
    I'll still need the input.mfa file, even though I won't be running BLAST
  --min_aa_coverage integer 1 to 100
  --min_aa_identity integer 1 to 100
  --min_aa_similarity integer 1 to 100
  -t tempdir (default is a subdir in /tmp)
  "
}
