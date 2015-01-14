#!/bin/sh

# download
wget --continue 'http://www.mediterranee-infection.com/arkotheque/client/ihumed/_depot_arko/articles/304/arg-annot-database_doc.zip'
if [ $? -gt 0 ]; then exit 1; fi;
unzip -o arg-annot-database_doc.zip && rm arg-annot-database_doc.zip
if [ $? -gt 0 ]; then exit 1; fi;

# format
legacy_blast.pl formatdb -p T -i Database_amino_acid_File.txt -t arg-annot -n arg-annot
if [ $? -gt 0 ]; then exit 1; fi;

