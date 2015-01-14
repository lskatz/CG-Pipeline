#!/usr/bin/env perl

# Parse uniprot XML file and output DB3 file with concise information that we need.
# This file is run during installation
# Author: Andrey Kislyuk
# Modifications: Lee Katz, Jay Humphrey

use strict;
use warnings;
use FindBin;
use lib "$FindBin::RealBin/../lib";
$ENV{PATH} = "$FindBin::RealBin:".$ENV{PATH};
use AKUtils qw(logmsg);
use File::Basename;

# debugging
use Data::Dumper;
#use Smart::Comments;

# database IO
use XML::LibXML::Reader;
use DBI; # SQLite
my ($dbfile, $evidence_dbfile) = ("cgpipeline.sqlite", "cgpipeline.evidence.sqlite");

# threads
use threads;
use threads::shared;
use Thread::Queue;
my $numCpus=AKUtils::getNumCPUs();
$numCpus=1; #debugging

# this gets you past a weird bug
# http://ubuntuforums.org/showthread.php?t=175050
##system("export MALLOC_CHECK_=0");
#if(!defined($ENV{MALLOC_CHECK_})){
  #die "Error: ENV variable MALLOC_CHECK_ must be 0\n";
#}
#$ENV{MALLOC_CHECK_}=0;

# globals
my $recordCount=0;
my $stick : shared; # YOU MAY ONLY WRITE TO THE DB IF YOU HAVE THE STICK!!!
my $dbh;

$0=fileparse($0);
exit(main());

sub main{
  die("Usage: $0 uniprot_sprot.xml uniprot_trembl.xml") if @ARGV < 1;
  my @infiles = @ARGV;
  for (@infiles) { die("File $_ not found") unless -f $_; }

  createDatabase();

  # set up the queue and threads that read from the queue
  my $queue=Thread::Queue->new;
  my @thr;
  my @dbh;
  for my $t(0 .. $numCpus-1){
    $dbh[$t]=connectToDatabase();
    push(@thr,threads->create(\&xmlReaderThread,$queue,$dbh[$t]));
  }

  # read the XML, send data to the queue for the threads
  INFILE:
  foreach my $infile (@infiles) {
    logmsg "Processing XML in file $infile with $numCpus CPUs";
    my $file_size = (stat $infile)[7];
    my $reader = new XML::LibXML::Reader(location => $infile) or die "cannot read $infile\n";
    my $i=0;
    # read the XML, paying close attention to the elements we want
    while ($reader->read) {
      # this is an element we want
      if ($reader->name eq 'entry' and $reader->nodeType != XML_READER_TYPE_END_ELEMENT) {
        $i++; $recordCount++; 
        # check in every 1000 records
        if($i % 1000 == 0){
          # status update
          my $status="$infile: ["
            .int(100*$reader->byteConsumed/$file_size)
            ."%] Processed $i records"
            ;

          # queue update
          if($i % 10000 == 0){
            $status.=" (jobs still pending: ".$queue->pending.")";
            # wait if the queue is really large
            sleep(1) while($queue->pending>10000);
          }
          logmsg($status);

          # debug
          if($i>2000){ logmsg "DEBUGGING, skipping rest of the file";next INFILE;}
        }
        # enqueue the data we want to look at; one thread will dequeue the data
        $queue->enqueue($reader->readOuterXml);
        $reader->next; # skip subtree
      }
    }
    NEXT_INFILE:
    logmsg "Processed $i records, done with file $infile";
    $reader->close; # free up any reader resources
  }
  # let the threads finish off the queue
  logmsg "Processed $recordCount records, done";
  while(my $numPending=$queue->pending){
    logmsg "$numPending more left..";
    sleep 1;
  }
  # the queue is finished; let the threads finish off what they have
  # This can be accomplished by grabbing the stick, indicating that no one else has the stick
  { lock($stick); } 

  # close off the database
  # TODO close db
  logmsg "All records are now in the DB!";

  # terminate the threads
  $queue->enqueue(undef) for(0..$numCpus-1); # send TERM signals
  while(my $t=shift @thr){
    if($t->is_joinable){
      $t->join;
    }
    else{push(@thr,$t);}
  }

  return 0;
}

# This is the subroutine of each thread.
# It waits for a queued item and then processes it
sub xmlReaderThread{
  my($queue,$dbh)=@_;
  while(my $data=$queue->dequeue){
    processUniprotXMLEntry($data,$dbh);
  }

  return 1;
}

# processes an XML entry
sub processUniprotXMLEntry($) {
	my ($xml_string,$dbh) = @_;
	my $parser = XML::LibXML->new();
	my $entry = $parser->parse_string($xml_string);

  # initialize all %info keys
	my %info;
  my @infoKeys=qw(accession name dataset proteinName proteinType geneType geneName dbRefId);

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

  # clear out any undef values
  $info{$_}||="" for @infoKeys;

  # prep the records for the database; put the keys in the right order
  my @line; # for %uniprot_h
  my @line_evidence; # for uniprot_evidence_h
  push(@line, $info{$_}) for @infoKeys;
  foreach my $ref (@db_refs) {
    my @line_ev; 
    push(@line_ev, $$ref{$_}) for qw(accession dbRefId dbRefType dbRefName);
    push(@line_evidence,\@line_ev);
  }
  # store the record in the database
  # but first lock the database to just this thread
  {
    lock($stick);
    #my $dbh=connectToDatabase();
    $_=$dbh->quote($_) for (@line);
    my $sql="INSERT INTO uniprot VALUES (".join(",",@line).")";
    print "$sql\n";
    $dbh->do($sql) or die("FAILED:\n$sql\n$!");

    foreach my $line (@line_evidence){
      $_=$dbh->quote($_) for (@$line);
      my $sql="INSERT INTO uniprot_evidence VALUES (".join(",",@$line).")";
      print "$sql\n";
      $dbh->do($sql);
    }
  }
  return 1;
}

sub createDatabase{
  unlink $dbfile; unlink $evidence_dbfile;
  
  my $dbh=connectToDatabase();
  
  $dbh->do("CREATE TABLE uniprot (accession,name,dataset,proteinName,proteinType,geneType,geneName,dbRefId)");
  $dbh->do("CREATE TABLE uniprot_evidence (accession,dbRefId,dbRefType,dbRefName)");

  return 1;
}
sub connectToDatabase{
  my $dbh;
  my $dbOptions={
    RaiseError=>1,
    sqlite_unicode=>1,
  };
  $dbh=DBI->connect("dbi:SQLite:dbname=$dbfile","","",$dbOptions);
  return $dbh;
}

