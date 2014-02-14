#!/usr/bin/env perl

# Author: Lee Katz <lkatz@cdc.gov>

# checks for all pipeline prerequisites
# Note: I found perl prereqs by doing the shell script
#  grep -h '^use\s*' *.pl|sed 's/qw.*$//'|perl -lane 's/use|;|\s//g;print'|grep -v 'use lib'|sort|uniq
# Note: for programs I used
#  grep -oh -i 'system("[a-z0-9_]\+' *.pl|perl -lane 's/system\(["]//;print' | sort|uniq
# TODO consider color text output

use strict;
use warnings;
use Getopt::Long;
use FindBin;
use lib "$FindBin::RealBin/../lib";
#$ENV{PATH} = "$FindBin::RealBin:".$ENV{PATH};
use AKUtils qw/logmsg/;
use File::Basename;
use Data::Dumper;

my $settings={
  appname=>'cgpipeline',

  # prereqs. Deprecated in favor of the hash below.
  perllibs=>[qw(AKUtils BerkeleyDB Bio::Assembly::IO Bio::Assembly::Scaffold Bio::Perl Bio::Seq Bio::SeqIO Bio::Seq::Quality Bio::Seq::RichSeq Bio::SeqUtils Bio::Tools::GFF Bio::Tools::Run::StandAloneBlast CGPBase CGPipeline::SQLiteDB CGPipelineUtils Data::Dumper Date::Format File::Basename File::Copy File::Path File::Spec File::Temp FindBin Getopt::Long GTTmhmm List::Util LWP::Simple Math::Round strict Thread::Queue threads threads::shared warnings XML::LibXML::Reader XML::Quote HTML::TableExtract IO::Scalar Mail::Send)],
  executables=>[qw(spades.py gam-create gam-merge addRun amos2ace AMOScmp cat cp gzip gunzip ln minimus2 mkdir newAssembly newMapping rm runProject setRef sfffile toAmos toAmos_new touch tmhmm signalp bam2fastq vcfutils.pl bcftools fastqqc velveth velvetg VelvetOptimiser.pl tRNAscan-SE gmsn.pl long-orfs extract build-icm glimmer3 rnammer hmmsearch iprscan)],
  environmentalVariables=>[qw(TMHMMDIR)],

  # presence/absence codes
  code_missing=>0,
  code_present=>1,
  code_usable=>2,
};

# each prerequisite with their descriptions.  Single or no characters in the values for no description.
my $prereqs={
  perllibs=>{
    pipeline=>{AKUtils=>1,BerkeleyDB=>1,'Bio::Perl'=>1,'Bio::Tools::Run::StandAloneBlast'=>1,CGPBase=>1,CGPipelineUtils=>1,'Data::Dumper'=>1,'Date::Format'=>1,'File::Basename'=>1,'File::Copy'=>1,'File::Path'=>1,'File::Spec'=>1,'File::Temp'=>1,FindBin=>1,'Getopt::Long'=>1,'List::Util'=>1,'LWP::Simple'=>1,'Math::Round'=>1, 'Thread::Queue'=>1, threads=>1, 'threads::shared'=>1, 'XML::LibXML::Reader'=>1, 'HTML::TableExtract'=>1},
    assembly=>{},
    prediction=>{GTTmhmm=>1},
    annotation=>{'XML::Quote'=>1,'Mail::Send'=>1},
  },
  executables=>{
    pipeline=>{cat=>1,cp=>1,gzip=>1, gunzip=>1, ln=>1,mkdir=>1,rm=>1, touch=>1, },
    assembly=>{'spades.py'=>1,'gam-create'=>1,'gam-merge'=>1,addRun=>1,amos2ace=>1,AMOScmp=>1,minimus2=>1,newAssembly=>1, newMapping=>1, runProject=>1, setRef=>1, sfffile=>1, toAmos=>1, toAmos_new=>1,bam2fastq=>1,'vcfutils.pl'=>1,bcftools=>1,'fastqqc'=>1, velveth=>1, velvetg=>1, 'VelvetOptimiser.pl'=>1},
    prediction=>{tmhmm=>1,signalp=>1,'tRNAscan-SE'=>1, 'gmsn.pl'=>1, 'long-orfs'=>1, extract=>1, 'build-icm'=>1, 'glimmer3'=>1, rnammer=>1, 'legacy_blast.pl'=>1},
    annotation=>{hmmsearch=>1, iprscan=>1, 'legacy_blast.pl'=>1},
  },
  environmentalVariables=>{
    pipeline=>{},
    assembly=>{},
    prediction=>{TMHMMDIR=>1},
    annotation=>{},
  },
};

exit(main());

sub main{
  $settings=AKUtils::loadConfig($settings);
  GetOptions($settings,qw(help verbose));
  die usage() if($$settings{help});

  my $envProblems=checkEnvVars($settings);
  my $exeProblems=checkExecutables($settings);
  my $libProblems=checkPerlLib($settings);

  logmsg "$envProblems problems found with environmental variables";
  logmsg "$libProblems problems found with libraries";
  logmsg "$exeProblems problems found with executables";

  my $totalProblems=$exeProblems+$libProblems;
  return 0 if(!$totalProblems);
  return (1+$totalProblems); # add one so that it doesn't look like a generic error
}

sub checkEnvVars{
  my($settings)=@_;
  my $problems=0;
  while(my($module,$environmentalVariable)=each(%{$$prereqs{environmentalVariables}})){
    while(my($var,$description)=each(%$environmentalVariable)){
      my $presence_code=is_envVar_present($var,$settings);
      $problems+=reportPresenceStatus($var,$presence_code,$settings);
    }
  }

  ## Any special env variables?
  # Check the CGPipeline path by seeing if this script is in the path
  local $0=basename $0;
  my $isInPath=AKUtils::fullPathToExec($0,{warn_on_error=>1});
  #die "=>$0:$isInPath<=";
  if(!$isInPath){
    $problems++;
    logmsg "CGP path not found. Need to add ".$FindBin::RealBin." to your PATH like so:";
    logmsg "  export PATH=".$FindBin::RealBin.":\$PATH";
  }
  return $problems;
}

sub checkPerlLib{
  my($settings)=@_;
  my $problems=0;
  while(my($module,$lib)=each(%{$$prereqs{perllibs}})){
    while(my($libname,$description)=each(%$lib)){
      my $presence_code=is_perlLib_present($libname,$settings);
      $problems+=reportPresenceStatus($libname,$presence_code,$settings);
    }
  }
  return $problems;
}

sub checkExecutables{
  my($settings)=@_;
  my $problems=0;
  while(my($module,$executable)=each(%{$$prereqs{executables}})){
    while(my($exec,$description)=each(%$executable)){
      my $presence_code=is_executable_present($exec,$settings);
      $problems+=reportPresenceStatus($exec,$presence_code,$settings);
    }
  }
  return $problems;
}

sub reportPresenceStatus{
  my($name,$presence_code,$settings)=@_;
  my $problems=0;
  if($presence_code == $$settings{code_missing}){
    logmsg "$name not found";
    $problems++;
  } elsif ($presence_code == $$settings{code_present}){
    logmsg "$name present, but not readable";
    $problems++;
  } elsif ($presence_code == $$settings{code_usable}){
    logmsg "$name is good!" if($$settings{verbose});
  }
  return $problems;
}

sub is_executable_present{
  my ($exe,$settings)=@_;
  my $code=$$settings{code_missing};
  for ("", split(/:/, $ENV{PATH})) {
    if(-e $_."/".$exe){
      next if(-d $_."/".$exe); # a directory is not "present but unusable" because it's not the actual exec
      $code=$$settings{code_present};
      if(-x $_."/".$exe){
        $code=$$settings{code_usable};
        return $code;
      }
    }
  }
  return $code;
}
sub is_perlLib_present{
  my($lib,$settings)=@_;

  my $relativeLibPath=$lib;
  $relativeLibPath=~s/::/\//g;
  $relativeLibPath.=".pm";
  
  # find the perl library
  my $code=$$settings{code_missing};
  for my $dir(@INC){
    my $absolutePath="$dir/$relativeLibPath";
    if(-e $absolutePath){
      $code=$$settings{code_present};
      if(-r "$dir/$relativeLibPath"){
        $code=$$settings{code_usable};
        return $code;
      }
    }
  }
  return $code;
}

sub is_envVar_present{
  my($env,$settings)=@_;
  my $code=$$settings{code_missing};
  if(defined $ENV{$env}){
    $code=$$settings{code_usable};
  }
  return $code;
}

sub usage{
  "Checks to see if all CG-Pipeline prerequisites are present on the system, in your environment.
  This script will have an error code of 1 + the number of problems found, or 0 for no problems.
  This script's error code will be 1 for general errors.
  Usage: $0
    -h for help
    -v for verbose
  "
}
