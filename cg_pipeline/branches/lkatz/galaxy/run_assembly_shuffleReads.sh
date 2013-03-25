#!/bin/sh

PE1=$1
PE2=$2
out=$3
TMP=/tmp

tmp1=$TMP/`basename $PE1`.fastq;
tmp2=$TMP/`basename $PE2`.fastq;

# gunzip or soft link each fastq file
if [ `echo $PE1|grep ".gz$"` ]; then 
  gunzip -c $PE1 > $tmp1
else
  ln -s $PE1 $tmp1
fi;
if [ `echo $PE2|grep ".gz$"` ]; then 
  gunzip -c $PE2 > $tmp2
else
  ln -s $PE2 $tmp2
fi;

# run the command and then remove the temporary files
run_assembly_shuffleReads.pl $tmp1 $tmp2 > $out
rm $tmp1 $tmp2
