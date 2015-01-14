#!/usr/bin/env perl

# Parse uniprot XML file and output DB3 file with concise information that we need.
# Author: Andrey Kislyuk
# Author: Lee Katz (lkatz@cdc.gov)

use strict;
use FindBin;
use lib "$FindBin::RealBin/../lib";
$ENV{PATH} = "$FindBin::RealBin:".$ENV{PATH};
use AKUtils qw(logmsg);
use Getopt::Long;
use File::Basename;

use XML::LibXML::Reader;
use BerkeleyDB;

use threads;
use threads::shared;
use Thread::Queue;

$0=fileparse($0);

my $settings={
  logfile=>"$0.log",
  numcpus=>AKUtils::getNumCPUs(),
};
my (%uniprot_h, %uniprot_evidence_h); # need to be global
my $recordsWritten :shared=0;
my $recordToStartFrom :shared=0;

$SIG{INT}=sub{closeDbs(); die @_;};
exit(main());

sub main{

  # figure out what the input files are
  my @infiles;
  for(@ARGV){
    last if(/^\-/);
    die("File $_ not found") unless -f $_;
    push(@infiles,$_);
  }
  if(!@infiles){
    @infiles=glob("*.xml") or die "No filenames supplied and no XML files found.\n".usage();

    # just to make them hit enter if they are sure of themselves
    logmsg "Query: use infiles: ".join(", ",@infiles)." \n? ctrl-C to exit, enter to continue? $0 -h for help.\n?";
    my $discard=<STDIN>;
  }
  # getting input files: done.

  GetOptions($settings,qw(force! help! logfile=s numExpectedEntries skipXml));
  die(usage()) if $$settings{help};

  reportNumEntries(\@infiles,$settings) if($$settings{numExpectedEntries});

  my $printQueue=Thread::Queue->new;
  my $printer=threads->new(\&txtDbPrinter,$printQueue,$settings);

  my $xmlToDbQueue=Thread::Queue->new;
  my @xmlToDb;
  $xmlToDb[$_]=threads->new(\&xmlToDbWorker,$xmlToDbQueue,$printQueue,$settings) for(0..$$settings{numcpus}-1);

  $recordToStartFrom=readLogfile($settings);

  foreach my $infile (@infiles) {
    last if($$settings{skipXml});
    logmsg "Processing XML in file $infile...";

    my $file_size = (stat $infile)[7];

    my $reader = new XML::LibXML::Reader(location => $infile) or die "cannot read $infile\n";	

    my $i;
    my @outerXml=();
    while ($reader->read) {
      if ($reader->name eq 'entry' and $reader->nodeType != XML_READER_TYPE_END_ELEMENT) {
        $i++;
        if($recordsWritten < $recordToStartFrom){ # useful for continuing
          $recordsWritten++;
          next;
        }
        if($i % 10000 == 0){
          logmsg("[".int(100*$reader->byteConsumed/$file_size)."%] Processed $i records...");
          # keep the queue down below 10k
          while($xmlToDbQueue->pending > 9999 || $printQueue->pending > 9999){
            logmsg "Slowing the analysis because the queue for printing to the database or the queue for analyzing the XML is quite large:". $printQueue->pending." ".$xmlToDbQueue->pending." respectively";
            sleep 1 while($xmlToDbQueue->pending>999 || $printQueue->pending>999);
          }
        }

        # build the queue in an array because enqueue is CPU intensive
        push(@outerXml,$reader->readOuterXml);
        if(@outerXml>9999){
          $xmlToDbQueue->enqueue(@outerXml);
          @outerXml=();
        }

        $reader->next; # skip subtree
      }
      #if($i>100000){
      #  logmsg "DEBUG";last;
      #}
    }
    $xmlToDbQueue->enqueue(@outerXml);
    logmsg "Processed $i records, done with file $infile";
  }
  logmsg "Processed $recordsWritten records, done";

  $xmlToDbQueue->enqueue(undef) for(0..$$settings{numcpus}-1);
  $_->join for(@xmlToDb);
  $printQueue->enqueue(undef);
  my($dbFile,$evidenceFile)=$printer->join;
  txtToDb3($dbFile,$evidenceFile,$settings);

  return 0;
}

sub xmlToDbWorker{
  my($Q,$printQueue,$settings)=@_;
  my @printBuffer;
  my $bufferLength=1000;
  my $i=0;
  while(defined(my $xmlStr=$Q->dequeue)){
    my($accession,$uniprot_h,$uniprot_evidence_h)=processUniprotXMLEntry($xmlStr);
    push(@printBuffer,[$accession,$uniprot_h,$uniprot_evidence_h]);
    $i++;
    if($i % $bufferLength==0){
      $printQueue->enqueue(@printBuffer);
      @printBuffer=();
    }
  }
  $printQueue->enqueue(@printBuffer) if(@printBuffer);
  return 1;
}

sub processUniprotXMLEntry($) {
	my ($xml_string) = @_;
	my $parser = XML::LibXML->new();
	my $entry = $parser->parse_string($xml_string);

	my %info;
	$info{accession} = $entry->getElementsByTagName('accession')->[0]->firstChild->nodeValue;

	$info{dataset} = $entry->getElementsByTagName('entry')->[0]->attributes->getNamedItem('dataset')->nodeValue;
	$info{name} = $entry->getElementsByTagName('name')->[0]->firstChild->nodeValue;
	# NB: the field names don't reflect the real meaning
	if (my $p = $entry->getElementsByTagName('protein')) {
		if ($p->[0]->childNodes and $p->[0]->childNodes->[1]->childNodes) {
			$info{proteinName} = $p->[0]->childNodes->[1]->childNodes->[1]->firstChild->nodeValue;
		}
	}
	if (my $pe_tag = $entry->getElementsByTagName('proteinExistence')) {
		$info{proteinType} = $pe_tag->[0]->attributes->getNamedItem('type')->nodeValue;
	}

	if (my @genes = $entry->getElementsByTagName('gene')) {
		my $name_tags = $genes[0]->getElementsByTagName('name');
		$info{geneName} = $name_tags->[0]->firstChild->nodeValue;
		$info{geneType} = $name_tags->[0]->attributes->getNamedItem('type')->nodeValue;
	} elsif (my @gene_location_list = $entry->getElementsByTagName('geneLocation')) {
		if (my $name_tags = $gene_location_list[0]->getElementsByTagName('name')) {
			if ($name_tags->[0]->firstChild) {
				$info{geneName} = $name_tags->[0]->firstChild->nodeValue;
				if ($info{geneType} = $name_tags->[0]->attributes->getNamedItem('type')) {
					$info{geneType} = $name_tags->[0]->attributes->getNamedItem('type')->nodeValue;
				}
			}
		}
	}

	my @db_refs;
	foreach my $db_ref ($entry->getElementsByTagName('dbReference')) {
		my $db_ref_type = $db_ref->attributes->getNamedItem('type')->nodeValue;
		my $db_ref_id = $db_ref->attributes->getNamedItem('id')->nodeValue;

		$info{dbRefId} = $db_ref_id if $db_ref_type eq 'GeneId';

		my @property_list = $db_ref->getElementsByTagName('property');
		next unless @property_list;
		my $property_type = $property_list[0]->attributes->getNamedItem('type')->nodeValue;
		next if $property_type ne 'entry name';
		
		my $db_ref_name = $property_list[0]->attributes->getNamedItem('value')->nodeValue;

		push(@db_refs, {accession => $info{accession},
						dbRefId => $db_ref_id,
						dbRefType => $db_ref_type,
						dbRefName => $db_ref_name});
	}

  my($uniprot_h,$uniprot_evidence_h)=("",[]);
	s/\|/\\|/g for values(%info);
	my @l; push(@l, $info{$_}) for qw(accession name dataset proteinName proteinType geneType geneName dbRefId);
  $uniprot_h=join('|',@l);

  foreach my $ref (@db_refs) {
		s/\|/\\|/g for values(%$ref);
		my @l; push(@l, $$ref{$_}) for qw(accession dbRefId dbRefType dbRefName);
    push(@$uniprot_evidence_h,join('|', @l));
    #$$uniprot_evidence_h{$info{accession}} = join('|', @l);
  }

  return($info{accession},$uniprot_h,$uniprot_evidence_h);
}

sub txtToDb3{
  my($dbfile,$evidence_dbfile,$settings)=@_;
  logmsg "$dbfile=>DB3";
  my $dbcounter=0;
  open(UNIPROT,"<",$dbfile) or die "Could not open $dbfile: $!";
  openDbs();
  #my (%uniprot_h, %uniprot_evidence_h); # need to be global
  while(<UNIPROT>){
    chomp;
    my $accession=(split(/\|/,$_))[0];
    $uniprot_h{$accession}=$_;
    logmsg "Finished with ". ($dbcounter/1000)."k db records" if(++$dbcounter%1000000==0);
  }
  close UNIPROT;
  my $evidencecounter=0;
  logmsg "$evidence_dbfile=>DB3";
  open(UNIPROT_EVIDENCE,"<",$evidence_dbfile) or die "Could not open $evidence_dbfile: $!";
  while(<UNIPROT_EVIDENCE>){
    chomp;
    my $accession=(split(/\|/,$_))[0];
    $uniprot_evidence_h{$accession}=$_;
    logmsg "Finished with ".($evidencecounter/1000)."k evidence records" if(++$evidencecounter%1000000==0);
  }
  close UNIPROT_EVIDENCE;
  return 1;
}
sub txtDbPrinter{
  my($Q,$settings)=@_;
  my($dbfile,$evidence_dbfile)=qw(uniprotdb.txt uniprotevidence.txt);
  unlink ($dbfile, $evidence_dbfile,$$settings{logfile}) if($$settings{force});
  #logmsg "DEBUG"; return ($dbfile,$evidence_dbfile);

  open(DB,">>",$dbfile) or die "Could not open $dbfile:$!";
  open(EVIDENCE,">>",$evidence_dbfile) or die "Could not open $evidence_dbfile:$!";
  while(defined(my $queueItem=$Q->dequeue)){
    my @cache=($queueItem);
    if($Q->pending>9999){
      push(@cache,$Q->dequeue(9990)); # leave room for a possible undef and then some
    }
    for (@cache){
      my($accession,$uniprot_h,$uniprot_evidence_h)=@$_;
      print DB $uniprot_h."\n";
      print EVIDENCE $_."\n" for(@$uniprot_evidence_h);
      $recordsWritten++;
    }
  }
  close DB; close EVIDENCE;
  return ($dbfile,$evidence_dbfile);
  return 1;
}

sub flushdb{
  $|++; $|--; # flush output too
  logmsg "flushing the database";
  closeDbs();
  openDbs();
  $|++; $|--; # flush output too
}

sub openDbs{
  my ($dbfile, $evidence_dbfile) = ("cgpipeline.db3", "cgpipeline.evidence.db3");
  tie(%uniprot_h, "BerkeleyDB::Hash", -Filename => $dbfile,  -Flags => DB_CREATE, -Property => DB_DUP)
    or die "Cannot open file $dbfile: $! $BerkeleyDB::Error\n";
  tie(%uniprot_evidence_h, "BerkeleyDB::Hash", -Filename => $evidence_dbfile, -Flags => DB_CREATE, -Property => DB_DUP)
    or die "Cannot open file $evidence_dbfile: $! $BerkeleyDB::Error\n";
  return 1;
}

sub closeDbs{
  logmsg "Closing the database at $recordsWritten records written.";
  untie %uniprot_h;
  untie %uniprot_evidence_h;
  logProgress($$settings{logfile},$recordsWritten,$settings);
  return 1;
}

sub reportNumEntries{
  my ($file)=@_;
  logmsg "Counting the number of expected entries in the XML files with grep. This could take a few minutes.";
  system("grep -cH '<entry' ".join(" ",@$file));
  warn "Could not find the number of entries: $!\n  Skipping.\n" if $?;
}

sub readLogfile{
  my ($settings)=@_;
  # if the logfile doesn't exist, give it a second try
  sleep  1 if(!-f $$settings{logfile});
  return 0 if(!-f $$settings{logfile});
  # read the log file and see where to continue from
  open(LOG,$$settings{logfile}) or die "Could not read from $$settings{logfile} because $!";
  my $line=<LOG>;
  chomp($line);
  close LOG;

  my($timestamp,$record)=split(/\t/,$line);
  $record=0 if($record<1); # just to be sure
  my($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime($timestamp);
  $year+=1900;
  logmsg "The last write was on $year-$mon-$mday with $record records. I am finding this record number now, but it might take a while if you have gone far. To start over, use the -f flag";
  if($$settings{force}){
    $record=0;
    logmsg "  However, you are forcing a restart and so I will start with record $record.";
  }
  
  return $record;
}
sub logProgress{
  my($logfile,$recordsWritten,$settings)=@_;
  $recordsWritten=$recordToStartFrom if($recordToStartFrom>$recordsWritten);
  open(LOG,">$logfile") or die "Could not open $logfile for writing because $!";
  print LOG join("\t",time(),$recordsWritten);
  print LOG "\n";
  close LOG;

  return 1;
}

sub usage{
  "
  Parses uniprot XML file and outputs DB3 file with concise information
  Usage: $0 uniprot_sprot.xml uniprot_trembl.xml
  If no xml files are supplied, then ./*.xml will be used
  -h for help (this text)
  -l logfile name (default: $$settings{logfile})
  -f to force
    logfile and database files will be deleted and you will start over
  -num to give an approximate count of entries before reading the XML files
  --skipXml to skip over reading the XML files (assume that you only need to read the text files)
  ";
}

# Need to flush database on exit
END{
  print "END of script detected. Flushing the databases.\n";
  closeDbs();
}

