#!/usr/bin/env perl

# Author: Lee Katz <lkatz@cdc.gov>

# checks for all pipeline prerequisites
# Note: I found perl prereqs by doing the shell script
#  grep -h '^use\s*' *.pl|sed 's/qw.*$//'|perl -lane 's/use|;|\s//g;print'|grep -v 'use lib'|sort|uniq
# Note: for programs I used
#  grep -oh -i 'system("[a-z0-9_]\+' *.pl|perl -lane 's/system\(["]//;print' | sort|uniq
# TODO consider color output

use strict;
use warnings;
use Getopt::Long;
use FindBin;
use lib "$FindBin::RealBin/../lib";
#$ENV{PATH} = "$FindBin::RealBin:".$ENV{PATH};
use AKUtils qw/logmsg/;
use File::Basename;

my $settings={
  appname=>'cgpipeline',
  perllibs=>[qw(AKUtils BerkeleyDB Bio::Assembly::IO Bio::Assembly::Scaffold Bio::Perl Bio::Seq Bio::SeqIO Bio::Seq::Quality Bio::Seq::RichSeq Bio::SeqUtils Bio::Tools::GFF Bio::Tools::Run::StandAloneBlast CGPBase CGPipeline::SQLiteDB CGPipelineUtils Data::Dumper Date::Format File::Basename File::Copy File::Path File::Spec File::Temp FindBin Getopt::Long GTTmhmm List::Util LWP::Simple Math::Round strict Thread::Queue threads threads::shared warnings XML::LibXML::Reader XML::Quote)],
  executables=>[qw(addRun amos2ace AMOScmp cat cp gzip gunzip ln minimus2 mkdir newAssembly newMapping rm runProject setRef sfffile toAmos toAmos_new touch tmhmm signalp bam2fastq vcfutils.pl bcftools fastqqc velveth velvetg VelvetOptimiser.pl tRNAscan-SE)],
  environmentalVariables=>[qw(TMHMMDIR)],

  # presence/absence codes
  code_missing=>0,
  code_present=>1,
  code_usable=>2,
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
  return (1+$totalProblems);
}

sub checkEnvVars{
  my($settings)=@_;
  my $problems=0;
  for my $var(@{$$settings{environmentalVariables}}){
    my $presence_code=is_envVar_present($var,$settings);
    $problems+=reportPresenceStatus($var,$presence_code,$settings);
  }
  # check the CGPipeline path by seeing if this script is in the path
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
  for my $lib(@{$$settings{perllibs}}){
    my $presence_code=is_perlLib_present($lib,$settings);
    $problems+=reportPresenceStatus($lib,$presence_code,$settings);
  }
  return $problems;
}

sub checkExecutables{
  my($settings)=@_;
  my $problems=0;
  for my $exe(@{$$settings{executables}}){
    my $presence_code=is_executable_present($exe,$settings);
    $problems+=reportPresenceStatus($exe,$presence_code,$settings);
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
