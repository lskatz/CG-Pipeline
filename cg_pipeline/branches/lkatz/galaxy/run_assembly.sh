#!/bin/sh
#run_assembly.sh $reads $expectedGenomeSize $out

READS=$1
SIZE=$2
OUT=$3

run_assembly $READS -e $SIZE -o $OUT

