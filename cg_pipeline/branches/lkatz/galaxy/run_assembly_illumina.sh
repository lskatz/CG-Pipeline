#!/bin/sh
# run_assembly_illumina.sh $reads $expectedGenomeSize $out
# Calls run_assembly for illumina reads only.
# By calling run_assembly, all read optimizations and cleanings are invoked

READS=$1
SIZE=$2
OUT=$3

run_assembly --illumina $READS -e $SIZE -o $OUT 2>&1

