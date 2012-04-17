#!/usr/bin/env perl

use strict;
use FindBin;
use lib "$FindBin::RealBin/../lib";
$ENV{PATH} = "$FindBin::RealBin:".$ENV{PATH};
use AKUtils qw(logmsg);
use Data::Dumper;
use File::Spec;
use Bio::Perl;
use Bio::SeqFeature::Generic;

use threads;
use Thread::Queue;

#use GTTmhmm;

# If TMHMM isn't set, it's the directory one level up from the bin directory
if(!$ENV{TMHMMDIR}){
  $ENV{TMHMMDIR} = `dirname \`which tmhmm\``;
  chomp($ENV{TMHMMDIR});
  $ENV{TMHMMDIR}=~s/\/\w+$//i;
  logmsg "Warning: ENV variable for TMHMMDIR is not set.  Setting it here as $ENV{TMHMMDIR}. Avoid this message in the future by adding the following line to ~/.bashrc\n  export TMHMMDIR=$ENV{TMHMMDIR}";
} else {
  logmsg "TMHMMDIR is set as $ENV{TMHMMDIR}";
}

exit(main());

sub main{
  my $settings={
    appname => 'cgpipeline',
  };
  die("Usage: $0 input.faa") unless -f $ARGV[0];
  $settings=AKUtils::loadConfig($settings);
  $$settings{numcpus}||=AKUtils::getNumCPUs();
  #logmsg "DEBUGGING numcpus=1"; $$settings{numcpus}=1;

  my $inseq = Bio::SeqIO->new(-file => "<$ARGV[0]", -format => 'fasta');

  # set up threads
  my $inQ=Thread::Queue->new;
  my $printOutQ=Thread::Queue->new;
  my $printStatsQ=Thread::Queue->new;
  my @thr;
  for(0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(\&tmhmmWorker,$inQ,$printOutQ,$printStatsQ,$settings);
  }
  my $printStatsThr=threads->new(\&printStatsWorker,$printStatsQ,$settings);
  my $printOutThr=threads->new(\&printOutWorker,$printOutQ,$settings);

  logmsg "Running TMHMM on proteins in $ARGV[0]...";
  my $i;
  while (my $seq = $inseq->next_seq) {
    $i++; logmsg "Enqueued $i proteins" if $i % 100 == 0;
    my $fastaEntry=seqToString($seq,"fasta",$settings);
    $inQ->enqueue($fastaEntry);
    if(0 && $i>30){
      logmsg "DEBUGGING quitting after 15 proteins";
      last;
    }
  }
  $inseq->close;
  logmsg "Done enqueuing $i proteins. Sending termination signal to threads.";
  $inQ->enqueue(undef) for(@thr);
  for(@thr){
    logmsg "Waiting on thread ".$_->tid;
    while($inQ->pending>$$settings{numcpus}+1){
      sleep 10;
      logmsg "Waiting to process ".$inQ->pending." sequences";
      # check the status every 10 seconds to shave some time off
      for(1..5){
        last if($inQ->pending<$$settings{numcpus}+1);
        sleep 10;
      }
    }
    $_->join;
  }
  logmsg "Done processing. Terminating print worker bees";
  $printStatsQ->enqueue(undef);
  $printOutQ->enqueue(undef);
  for($printStatsThr,$printOutThr){
    $_->join;
  }
  logmsg "Processed $i proteins, done";

  return 0;
}

sub tmhmmWorker{
  my($inQ,$printOutQ,$printStatsQ,$settings)=@_;
  my $tid=threads->tid;
  my $tempdir=File::Spec->tmpdir()."/TMHMM.$tid.$$";

  # The TMHMM runner and parser had to be modified to allow parsing of the "inside"/"outside" markers and stats.
  while(defined(my $fastaEntry=$inQ->dequeue)){
    # convert fastaEntry back to a SeqObj
    my $stringfh;
    open($stringfh,"<",\$fastaEntry) or die "Could not open string for reading because $!";
    my $seq=Bio::SeqIO->new(-fh=>$stringfh,-format=>"fasta")->next_seq;

    my ($feats, $stats) = runTmhmm($fastaEntry);
    next if @$feats < 2; # a transmembrane motif results in at least 3 features (inside-transmembrane-outside)

    foreach my $feat (@$feats) {
      my $type = $feat->primary_tag;
      $type =~ s/TMhelix/Transmembrane Helix/;
      $printOutQ->enqueue(join("|", ($seq->id, $type, $feat->start, $feat->end))."\n");
    }

    my $tmhmm_measures = ["Length: ", "Number of predicted TMHs: ", "Exp number of AAs in TMHs: ", "Exp number, first 60 AAs: ", "Total prob of N-in: "];
    my @l = ( $seq->id );
    
    for (@$tmhmm_measures) { push(@l, $$stats{$_}); }
    $printStatsQ->enqueue(join("|", @l)."\n");
  }

  return 1;
}

sub printStatsWorker{
  my($printQ,$settings)=@_;
  open(TMHMM_STATS_SQL, '>', "$ARGV[0].tmhmm.sql") or die "Could not open file for writing: $!";
  while(defined(my $line=$printQ->dequeue)){
    #logmsg "Printing: $line";
    print TMHMM_STATS_SQL $line;
  }
  close TMHMM_STATS_SQL;
  return 1;
}
sub printOutWorker{
  my($printQ,$settings)=@_;
  open(TMHMM_OUT_SQL, '>', "$ARGV[0].tmhmm_location.sql") or die "Could not open file for writing: $!";
  while(defined(my $line=$printQ->dequeue)){
    #logmsg "Printing: $line";
    print TMHMM_OUT_SQL $line;
  }
  close TMHMM_OUT_SQL;
  return 1;
}

sub seqToString{
  my($seq,$format,$settings)=@_;
  $format||="fasta";
  my $str;
  my $strfh;
  open($strfh,">",\$str) or die "Could not open string for writing because $!";
  my $seqout=Bio::SeqIO->new(-format=>$format,-fh=>$strfh);
  $seqout->write_seq($seq);
  close FASTAOUT;
  return $str;
}

# mostly taken from bioperl
# parameter: string of a fasta entry
sub runTmhmm{
  my($fasta,$settings)=@_;

  my $program_dir=$ENV{TMHMMDIR} || '';
  my $tempdir=AKUtils::mktempdir($settings);
  my $infile="$tempdir/in.fasta";
  my $outfile="$tempdir/out.tmhmm";
  open(FASTA,">$infile") or die "Could not open $infile for writing: $!";
  print FASTA $fasta;
  close FASTA;

  # run the program
  my $str=AKUtils::fullPathToExec("tmhmm");
  $str.=" --noplot -basedir=$program_dir -workdir=$tempdir $infile > $outfile";
  system("$str");
  die "Error with TMHMM: $!" if $?;

  # read the results
  my @tmhmm_feat; my %stats;
  open(TMHMMOUT,$outfile) or die "Could not read the results file: $!";
  while(my $line=<TMHMMOUT>){
    if ( $line =~ m/^(\S+)\s+(\S+)\s+(\w+)\s+(\d+)\s+(\d+)$/i ) {
      push(@tmhmm_feat, Bio::SeqFeature::Generic->new(
        -primary => $3,
        -seq_id  => $1,
        -source  => $2,
        -start   => $4,
        -end     => $5,
      ));
    } elsif ($line =~ /^#/) {
      my $tmhmm_measures = ["Length: ", "Number of predicted TMHs: ", "Exp number of AAs in TMHs: ", "Exp number, first 60 AAs: ", "Total prob of N-in: "];
      foreach my $measure (@$tmhmm_measures) {
        if ($line =~ /$measure\s*(.+)/) {
          $stats{$measure} = $1;
        }
      }
    }
  }
  close TMHMMOUT;
  return(\@tmhmm_feat,\%stats);
}

