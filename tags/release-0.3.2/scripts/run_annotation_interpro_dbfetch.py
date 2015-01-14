#! /usr/bin/python

# Author: Matthew Hagen

import urllib
from xml.dom import minidom
import os
import sys
import re
from optparse import OptionParser

#parse input arguments
usage = "usage: %prog [option] arg"
parser = OptionParser(usage=usage)
parser.add_option("-f", action="store", type="string", dest="filename", help="should be teh file that outputs interpro_domain_unique.sql")
parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True, help="don't print status messages to stdout")

(options, args) = parser.parse_args()

if len(args) != 0:
	sys.exit("Error parsing arguments, "+usage);
	
phandle = open(options.filename)

interpro_url = "http://www.ebi.ac.uk/cgi-bin/dbfetch?db=interpro&id=%s&format=interproxml"

db_filename 				= "interpro_db2.sql"
db_link_filename 			= "interpro_db_link2.sql"
pubmed_link_filename 			= "interpro_pubmed_link2.sql"
go_link_filename 			= "interpro_go_link2.sql"
go_filename 				= "interpro_go2.sql"
domain_filename 			= "interpro_domain_nen_nen_new.sql"
error_filename                          = "error_new.sql"

db_filename_handle 			= open( db_filename, 'w' )
db_link_filename_handle 		= open( db_link_filename, 'w' )
pubmed_link_filename_handle 		= open( pubmed_link_filename, 'w' )
go_link_filename_handle 		= open( go_link_filename, 'w' )
go_filename_handle 			= open( go_filename, 'w' )
domain_filename_handle 			= open( domain_filename, 'w' )
error_filename_handle                   = open(error_filename, 'w')

#parse input paths to the sequences 
for line in phandle:				
	ipr_id = line.rstrip()
	sys.stderr.write("Processing id "+ipr_id+"\n")
	
	url = interpro_url % ipr_id

	try:
		minidom.parse(urllib.urlopen(url))
		xmldoc = minidom.parse(urllib.urlopen(url))
	except:
		error_filename_handle.write(ipr_id+" "+"not found" +"\n")

	# parse each of the IDs received from interproscan to the database [interpro_url]
	# so as to use the xml data about each sequence to retrieve information about that ID. 

	interproList = xmldoc.getElementsByTagName('interpro')

	interpro = interproList[0]
	
	id = interpro.attributes["id"].value.encode('iso-8859-1')

	type = interpro.attributes["type"].value.encode('iso-8859-1')

	name = interpro.getElementsByTagName('name')

	name = name[0].firstChild.data.encode('iso-8859-1')
	abstract = interpro.getElementsByTagName('abstract')[0]
	
	abstractText = ""
	abstractList = abstract.childNodes
	for abstract in abstractList:
		if (abstract.nodeType == abstract.TEXT_NODE):
			abstractText += abstract.data.encode('iso-8859-1')

	abstractText = re.sub("\s+", " ", abstractText).strip()	

	output = id + "|" + name + "|" + type + "|" + abstractText
	domain_filename_handle.write(output + "\n")

	classList = interpro.getElementsByTagName('class_list')
	
	if (len(classList) > 0):
	
		classList = classList[0].getElementsByTagName('classification')
	
		for classification in classList:
	
			classId = classification.attributes["id"].value
			classType = classification.attributes["class_type"].value
			
			category = classification.getElementsByTagName('category')[0].firstChild.data
			description = classification.getElementsByTagName('description')[0].firstChild.data
			
			if (classType == "GO"):
				output = classId + "|" + category + "|" + description
				go_filename_handle.write(output + "\n")
				output = id + "|" + classId
				go_link_filename_handle.write(output + "\n")

	pubList = interpro.getElementsByTagName('pub_list')
	
	if (len(pubList) > 0):	
		pubList = pubList[0].getElementsByTagName('publication')
		
		for publication in pubList:
			pubmedid = publication.getElementsByTagName('db_xref')
			if (len(pubmedid) == 0):
				continue
				
			pubmedid = pubmedid[0].attributes['dbkey'].value
			output = id + "|" + pubmedid
			
			pubmed_link_filename_handle.write(output + "\n")
	
	memberList = interpro.getElementsByTagName('member_list')
		
	if (len(memberList) > 0):
		memberList = memberList[0].getElementsByTagName('db_xref')
		
		for member in memberList:
			db = member.attributes['db'].value
			dbKey = member.attributes['dbkey'].value
			name = member.attributes['name'].value
		
			output = id + "|" + dbKey
			db_link_filename_handle.write(output + "\n")
			output = dbKey + "|" + name + "|" + db
			db_filename_handle.write(output + "\n")
	
db_filename_handle.close()	
db_link_filename_handle.close()
pubmed_link_filename_handle.close()
go_link_filename_handle.close()
go_filename_handle.close()
domain_filename_handle.close()
error_filename_handle.close()
