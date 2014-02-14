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

  # presence/absence codes
  code_missing=>0,
  code_present=>1,
  code_usable=>2,
};

# each prerequisite with their descriptions.  Single or no characters in the values for no description.
my $prereqs={
  perllibs=>{
    pipeline=>{AKUtils=>"CGP Module",BerkeleyDB=>1,'Bio::Perl'=>"BioPerl",'Bio::Tools::Run::StandAloneBlast'=>1,CGPBase=>"CGP Module",CGPipelineUtils=>"CGP Module",'Data::Dumper'=>1,'Date::Format'=>1,'File::Basename'=>1,'File::Copy'=>1,'File::Path'=>1,'File::Spec'=>1,'File::Temp'=>1,FindBin=>1,'Getopt::Long'=>"For parsing options",'List::Util'=>1,'LWP::Simple'=>1,'Math::Round'=>1, 'Thread::Queue'=>1, threads=>1, 'threads::shared'=>1, 'XML::LibXML::Reader'=>1, 'HTML::TableExtract'=>1},
    assembly=>{},
    prediction=>{GTTmhmm=>"TMHMM module"},
    annotation=>{'XML::Quote'=>"Required by InterProScan",'Mail::Send'=>"Required by InterProScan"},
  },
  executables=>{
    pipeline=>{cat=>1,cp=>1,gzip=>1, gunzip=>1, ln=>1,mkdir=>1,rm=>1, touch=>1, },
    assembly=>{'spades.py'=>"SPAdes assembler",'gam-create'=>"NGS-GAM assembly merger",'gam-merge'=>"NGS-GAM assembly merger",addRun=>"Newbler",amos2ace=>"AMOS",AMOScmp=>"AMOS",minimus2=>"AMOS",newAssembly=>"Newbler", newMapping=>"Newbler", runProject=>"Newbler", setRef=>"Newbler", sfffile=>"Newbler", toAmos=>"AMOS", toAmos_new=>"AMOS",bam2fastq=>"Exporting fastq from bam files; found in cg_pipeline/etc",'vcfutils.pl'=>"Samtools",bcftools=>"Samtools",'fastqqc'=>"AMOS", velveth=>"Velvet", velvetg=>"Velvet", 'VelvetOptimiser.pl'=>"Velvet"},
    prediction=>{tmhmm=>"TMHMM transmembrane helix predictor",signalp=>"Signal Peptides",'tRNAscan-SE'=>"tRNA prediction", 'gmsn.pl'=>"GeneMark", 'long-orfs'=>"Glimmer3", extract=>"Glimmer3", 'build-icm'=>"Glimmer3", 'glimmer3'=>"Glimmer3", rnammer=>"RNA prediction", 'legacy_blast.pl'=>"BLAST+", blastp=>"BLAST+"},
    annotation=>{hmmsearch=>1, iprscan=>"InterProScan", 'legacy_blast.pl'=>"BLAST+", blastp=>"BLAST+"},
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
      $problems+=reportPresenceStatus($var,$description,$module,$presence_code,$settings);
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
      $problems+=reportPresenceStatus($libname,$description,$module,$presence_code,$settings);
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
      $problems+=reportPresenceStatus($exec,$description,$module,$presence_code,$settings);
    }
  }
  return $problems;
}

sub reportPresenceStatus{
  my($name,$desc,$module,$presence_code,$settings)=@_;
  $module=uc($module);

  # set up the description of the thing to report presence/absence
  my $description="no description";
  if($desc && length($desc)>1){ # must have >1 character to be informative
    $description=$desc;
  }

  my $problems=0;
  if($presence_code == $$settings{code_missing}){
    logmsg "$name ($module - $description) not found";
    $problems++;
  } elsif ($presence_code == $$settings{code_present}){
    logmsg "$name ($module - $description) present, but not readable";
    $problems++;
  } elsif ($presence_code == $$settings{code_usable}){
    logmsg "$name ($module - $description) is good!" if($$settings{verbose});
  }
  return $problems;
}

sub is_executable_present{
  my ($exe,$settings)=@_;
  my $code=$$settings{code_missing};
  for ("", split(/:/, $ENV{PATH})) {
    my $possibleExec=$_."/".$exe;
    # If the possibleExec exists and is a file...
    if(-e $possibleExec && -f $possibleExec){
      $code=$$settings{code_present};
      if(-x $possibleExec){
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
  EXAMPLE
    $0 | grep ASSEMBLY # see what is missing to perform assembly
  "
}
