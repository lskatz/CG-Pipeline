#!/bin/sh
# run_assembly_illumina.sh $reads $expectedGenomeSize $out
# Calls run_assembly for illumina reads only.
# By calling run_assembly, all read optimizations and cleanings are invoked

READS=$1
SIZE=$2
OUT=$3

. `dirname $0`/cgpipelineGalaxyrc

run_assembly --illumina $READS -e $SIZE -o $OUT 2>&1
if [ $? -gt 0 ]; then
  echo "ERROR with Illumina assembly. See the log for details." 1>&2
  exit 1
fi

