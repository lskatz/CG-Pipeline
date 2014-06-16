#!/usr/bin/env perl

# CGPipelineUtils: more subroutines

package CGPipelineUtils;
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
use AKUtils qw/alnum printSeqsToFile/;
use Data::Dumper;

local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };

our @ISA = "Exporter";
our @methods = qw();
our %EXPORT_TAGS = (all => [@methods]);
Exporter::export_ok_tags('all');

sub logmsg {my $FH = $AKUtils::LOG || *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

# Run and parse SignalP v4 output.
# Both the runner and the parser in BioPerl were found inadequate.
# Usage: $predictions = getSignalPPredictions($seqs, {signalp_org_type => 'gram-', ...})
#     $seqs = {seq1name => seq1, ...}
#     $predictions = {seq1name => [ {type=>..., ...}, ...], ...}
sub getSignalPPredictionsV4($$) {
  my ($seqs, $settings) = @_;
  die("No input sequences supplied") if ref($seqs) ne 'HASH';
  die("Setting signalp_org_type must be set to gram-, gram+ or euk") if $$settings{signalp_org_type} !~ /^(gram\-|gram\+|euk)$/;
  die("Invalid value of setting signalp_pred_method, must be one of nn (neural networks), hmm (hidden Markov models), nn+hmm")
    if defined $$settings{signalp_pred_method} and $$settings{signalp_pred_method} !~ /^(nn|hmm|nn\+hmm)$/;
  $$settings{signalp_trunc_length} = 70 if not defined $$settings{signalp_trunc_length};
  $$settings{numcpus}||=AKUtils::getNumCPUs();
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
  $signalp_opts .= " -c $$settings{signalp_trunc_length}";

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
    # try one more time just in case
    if($?){
      logmsg "SignalP ERROR: $!";
      logmsg "Trying one more time.";
      sleep 1; # pause just in case something is conflicting with an input temporary file or whatever
      system($invoke_string);
      die("Error running signalp: $!") if $?;
    }

    my $result = loadSignalPPredictions($signalpOutfile);
    $predictions = {%$predictions, %$result};
  }
  return $predictions;
}

sub loadSignalPPredictions($) {
  my ($predfile) = @_;
  my @predictions;
  my $pred_text;
  open(PRED, '<', $predfile) or die("Could not open file $predfile for reading: ".$!);
  { local $/; $pred_text = <PRED>; }
  my @pred_blocks=grep(!/^\s*$/,split(/^#/,$pred_text));
  @pred_blocks=map("#$_",@pred_blocks); # put back on the hash sign for parsing
  close PRED;

  my %prediction;
  for my $pred_block (@pred_blocks){
    my $name;
    for my $line(split /\n/,$pred_block){
      my $prediction={type=>"nn"};
      if($line=~/(max\. ([CSY])|mean S| D )/){
        # split on 2+ spaces because some columns have a space
        my @l=grep(!/^\s*$/,split(/\s{2,}/, $line));
        for (qw(measure position value cutoff decision)){
          $$prediction{$_}=shift @l;
        }
      }
      if($line=~/Name=(\S+)/){
        $name=$1;
      }
      ($$prediction{start}, $$prediction{end}) = split(/-/, $$prediction{position});
      next if(!$$prediction{measure});
      push(@predictions,$prediction);
    }
    if($name){
      $prediction{$name}=\@predictions;
    }
  }
  return \%prediction;
}

sub getProdigalPredictions{
  my($assembly,$settings)=@_;

  # run Prodigal
  my $asm=join(" ",@$assembly);
  my $gff=`run_prediction_prodigal.pl $asm -t $$settings{tempdir}`; die if $?;
  my @gff=split("\n",$gff);
  
  # read the gff file
  my %predictions;
  for(@gff){
    next if(/^\s*#/);
    my($seqname,$predictor,$type,$lo,$hi,$score,$strand,undef,$attribute)=split /\t/;
    my $p={lo=>$lo,hi=>$hi,predictor=>$predictor,seqname=>$seqname,strand=>$strand,type=>$type};
    if($strand eq '+'){
      $$p{start}=$lo;
      $$p{stop}=$hi;
    } else {
      $$p{start}=$hi;
      $$p{stop}=$lo;
    }
    push(@{$predictions{$seqname}},$p);
  }
  return \%predictions;
}

1;

