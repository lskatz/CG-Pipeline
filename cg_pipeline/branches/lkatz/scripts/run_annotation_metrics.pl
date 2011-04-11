#!/usr/bin/env perl

# run_annotation_metrics.pl: finds metrics of a given annotation genbank file
# Author: Lee Katz (lskatz@gatech.edu)

package PipelineRunner;
my ($VERSION) = ('$Id: $' =~ /,v\s+(\d+\S+)/o);

my $settings = {
    appname => 'cgpipeline',
};
my $stats;

use strict;
no strict "refs";
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
use Data::Dumper;
use Bio::Perl;

$0 = fileparse($0);
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {my $FH = $FSFind::LOG || *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}

exit(main());

sub main() {
  $settings = AKUtils::loadConfig($settings);
  die("  Usage: $0 -p projectName\n") if @ARGV<1;

  my @cmd_options=('output=s','project=s');
  GetOptions($settings, @cmd_options) or die;

  my $project=$$settings{project};
  my $project=File::Spec->rel2abs($project);
  die("Project $project not found") unless -d $project;

  my $metrics=annotationMetrics($project,$settings);
  # TODO print metrics in an understandable and parsable way
  print join("\t","Project",$project)."\n";
  foreach my $metric (keys %$metrics){
    print join("\t",$metric,$$metrics{$metric},($$metrics{$metric}/$$metrics{CDS}*100))."\n";
  }
  return 0;
}

sub annotationMetrics($$){
  my($project,$settings)=@_;
  my $metrics; # hash of feature metrics
  my $metricsSub; # hash of tag metrics
  my %thingsToCount=("vfdb_id"=>1);
  my $count;
  
  my %seenCDS;
  my $file="$project/annotation.gb";
  my $seqio=Bio::SeqIO->new(-file=>$file,-format=>'genbank');
  while(my $seq=$seqio->next_seq){
    FEATURE: for my $feat ($seq->get_SeqFeatures){
      $$metrics{$feat->primary_tag}++;
      #next if($seenFeature{$feat->primary_tag});
      #$seenFeature{$feat->primary_tag}++;
      for my $tag ($feat->get_all_tags){
        $$metricsSub{$tag}++;
        # count certain values
        if($thingsToCount{$tag}){
          $$metrics{$tag}++;
        }
        for my $value($feat->get_tag_values($tag)){
          # count certain misc features
          if($feat->primary_tag=~/^(misc_structure|misc_feature)$/){
            if($value=~/Transmembrane Helix/i){
              $$metrics{transmembraneHelix}++;
              next FEATURE;
            }
          }

          # see which genes have a name or not
          if($feat->primary_tag eq 'gene'){
            if($tag eq 'gene' && $value ne ''){
              $$metrics{characterized}++;
              next FEATURE;
            }
          }
          # CDS stats
          if($feat->primary_tag eq 'CDS'){
            # how many are "putative" or "hypothetical"
            if($tag eq 'product' && $value=~/(putative|hypothetical)/){
              $$metrics{putativeCDS}++ if($1 eq 'putative');
              $$metrics{hypotheticalCDS}++ if($1 eq 'hypothetical');
              next FEATURE;
            }
          }
        }
      }
    }
  }
  # find how many genes have at least some annotation
  my $totalGenes;
  my $genesAnnotated=0;
  my $foundAProteinProduct=0;
  $seqio=Bio::SeqIO->new(-file=>$file,-format=>'genbank');
  while(my $seq=$seqio->next_seq){
    FEATURE: for my $feat ($seq->get_SeqFeatures){
      if($feat->primary_tag eq 'gene'){
        $totalGenes++;
        $foundAProteinProduct=0;
      }
      for my $tag ($feat->get_all_tags){
        if(!$foundAProteinProduct && $tag eq 'product'){
          $genesAnnotated++;
          $foundAProteinProduct=1;
        }
      }
    }
  }
  $$metrics{genesAnnotated}=$genesAnnotated;
  #$metrics{raw}=getFileLineNumbers($project); # only if something isn't captured right--this will be a hack
  #print Dumper $metricsSub;print Dumper $metrics;exit;

  return $metrics;
}

# count the number of genes in all the sql files
sub getFileLineNumbers{
  my ($project)=@_;
  my $lineNumber;
  my @file=glob("$project/annotation/*.sql");
  for my $f (@file){
    my($filename)=fileparse($f);
    $$lineNumber{$filename}=`sed 's/|/\t/g' $f|cut -f 1|sort|uniq |wc -l`;
  }
  chomp(%$lineNumber);
  return $lineNumber;
}

