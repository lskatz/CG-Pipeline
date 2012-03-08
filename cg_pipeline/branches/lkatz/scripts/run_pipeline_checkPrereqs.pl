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
$ENV{PATH} = "$FindBin::RealBin:".$ENV{PATH};
use AKUtils qw/logmsg/;

my $settings={
  appname=>'cgpipeline',
  perllibs=>[qw(AKUtils BerkeleyDB Bio::Assembly::IO Bio::Assembly::Scaffold Bio::Perl Bio::Seq Bio::SeqIO Bio::Seq::Quality Bio::Seq::RichSeq Bio::SeqUtils Bio::Tools::GFF Bio::Tools::Run::StandAloneBlast CGPBase CGPipeline::SQLiteDB CGPipelineUtils Data::Dumper Date::Format File::Basename File::Copy File::Path File::Spec File::Temp FindBin Getopt::Long GTTmhmm List::Util LWP::Simple Math::Round strict Thread::Queue threads threads::shared warnings XML::LibXML::Reader)],
  executables=>[qw(addRun amos2ace AMOScmp cat cp gunzip ln minimus2 mkdir nesoni newAssembly newMapping rm runProject setRef sfffile toAmos toAmos_new touch)],

  # presence/absence codes
  code_missing=>0,
  code_present=>1,
  code_usable=>2,
};

exit(main());

sub main{
  $settings=AKUtils::loadConfig($settings);
  GetOptions($settings,qw(help errorsOnly));

  my $exeProblems=checkExecutables($settings);
  my $libProblems=checkPerlLib($settings);

  logmsg "$libProblems problems found with libraries";
  logmsg "$exeProblems problems found with executables";
  return 0;
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
    logmsg "$name is good!" if(!$$settings{errorsOnly});
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

sub usage{
  "Checks to see if all CG-Pipeline prerequisites are present on the system, in your environment.
  Usage: $0
    -h for help
    -e for print only errors
  "
}
