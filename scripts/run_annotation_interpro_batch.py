#!/usr/bin/env python

# Author: Matthew Hagen

# Import WSDL package
import time
from SOAPpy import WSDL
from Bio import Fasta 
from xml.dom import minidom

from Bio.Alphabet import IUPAC
from Bio import Translate
from Bio import SeqIO
import os
import sys
import re
from optparse import OptionParser

# Process command-line options
usage = "usage: %prog [option] arg"
parser = OptionParser(usage=usage)
#parser.add_option("-f", action="store", type="string", dest="infile", help="Reads name of the file storing absolute paths for gene prediction results")
parser.add_option("--interpro_evidence_sql_outfile", action="store", type="string",
	dest="interpro_evidence_sql_outfile", help="Reads name of the file storing absolute paths for gene prediction results")
parser.add_option("--interpro_domain_sql_outfile", action="store", type="string",
	dest="interpro_domain_sql_outfile", help="Reads name of the file storing absolute paths for gene prediction results")
parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True, help="don't print status messages to stdout")
parser.add_option('--reporting_email', action="store", type="string", dest="reporting_email", help='E-mail address')

parser.add_option('--app', help='List of methods to run')
parser.add_option('--crc', action="store_true", help='Use IprMatches lookup')
parser.add_option('--seqtype', default='P', help='Sequence type')
#parser.add_option('--trlen', type=int, help='Minimum ORF length')
#parser.add_option('--trtable', type=int, help='Translation table to use')
parser.add_option('--goterms', action="store_true", help='Get GO terms')
parser.add_option('-f', '--outformat', help='Output format')
parser.add_option('-a', '--async', action="store_true", help='Async submission')
parser.add_option('-o', '--outfile', help='File name for results')
parser.add_option('--polljob', action="store_true", help='Get job result')
parser.add_option('--status', action="store_true", help='Get job status')
parser.add_option('-j', '--jobid', help='Job Id')
parser.add_option('--trace', action="store_true", help='Show SOAP messages')
parser.add_option('--sequence', help='Input sequence file name')
#parser.add_option('--handle_blast', help= 'handle blast outputs')
(options, args) = parser.parse_args()

if options.interpro_domain_sql_outfile is None: sys.exit("Argument interpro_domain_sql_outfile is required. "+usage)
if options.interpro_evidence_sql_outfile is None: sys.exit("Argument interpro_evidence_sql_outfile is required. "+usage)
if options.reporting_email is None: sys.exit("Argument reporting_email is required. "+usage)

if len(args) != 1: sys.exit("Input filename not supplied. "+usage)

infile = args[0]
if not os.path.exists(infile): sys.exit("Input file "+infile+" not found")

if options.outfile is None: options.outfile = infile+".interpro_batch.out"

# Create service interface
# iprscan WSDL documentation: http://www.ebi.ac.uk/Tools/webservices/services/interproscan
wsdlUrl = 'http://www.ebi.ac.uk/Tools/webservices/wsdl/WSInterProScan.wsdl'
server = WSDL.Proxy(wsdlUrl)

iprscan_params = {
	'email': options.reporting_email,
	'app': options.app,
	'seqtype': options.seqtype,
}

iprscan_params['async'] = 1
iprscan_params['goterms'] = 1
iprscan_params['crc'] = 1

wh_evidence = open(options.interpro_evidence_sql_outfile, 'w')
wh_domain = open(options.interpro_domain_sql_outfile, 'w')

# parse results to get piped outputs in sql format
def parseResult(result, handle_domain, handle_evidence, protein):
	xmldoc = minidom.parseString(result)
	
	interproList = xmldoc.getElementsByTagName('interpro')
	for interpro in interproList:
		id = interpro.attributes["id"].value

		if ( id != "noIPR"):
			handle_domain.write(id + "\n")
		matchList = interpro.getElementsByTagName('match')
		for match in matchList:
			id = match.attributes["id"].value

			dbname = match.attributes["dbname"].value
			locationList = match.getElementsByTagName('location')
			location = locationList[0]
			start = location.attributes["start"].value
			end = location.attributes["end"].value
			score = location.attributes["score"].value
			status = location.attributes["status"].value
			evidence = location.attributes["evidence"].value
			output = protein + "|" + id + "|" + dbname + "|" + start + "|" + end + "|" + score + "|" + status + "|" + evidence + "\n"
			handle_evidence.write(output)

#check if jobs are completed

def checkJobs(jobIds, wh_domain, wh_evidence):
	counter = 0
	for jobIdStr in jobIds:
		elements = jobIdStr.split(":")
		title = elements[0]
		jobId = elements[1]
		
		result = server.checkStatus(jobId)
	
		if (result != 'RUNNING') and (result != 'PENDING'):
			jobIds = jobIds[:counter] + jobIds[(counter+1):]
		
			sys.stderr.write("Job finished: "+result+": " + title + " | " + jobId + "\n")
		
			resultTypes = server.getResults(jobId)
			
			for resultType in resultTypes:
				if (resultType.type == "toolxml"):			
					filename = "interpro_" + title + ".xml"
					result = server.poll(jobId, resultType.type)
					parseResult(result, wh_domain, wh_evidence, title)
		else:
			counter += 1

	return jobIds

	
jobIds = []
packageFull = False

RH = open(infile, "r")
WH = open(options.outfile, "w")

for prot_record in Fasta.Iterator(RH, Fasta.RecordParser()):
	seqData = ">" + prot_record.title + "\n" + prot_record.sequence
	content = [{'type':'sequence', 'content':seqData}]
	packageFull = True

	while (packageFull):
		jobIds = checkJobs(jobIds, wh_domain, wh_evidence)

		if (len(jobIds) < 25):
			packageFull = False
						
		if (packageFull):
			time.sleep(15)

	sys.stderr.write("Submitted protein "+prot_record.title+" to InterProScan...\n")
	
	WH.write(">" + prot_record.title + "\n")

	jobId = server.runInterProScan(params=iprscan_params, content=content)
		
	jobIdStr = prot_record.title + ":" + jobId
	jobIds += [jobIdStr]
	WH.write(jobId + "\n")
RH.close()
WH.close()

wh_evidence.close()
wh_domain.close()
