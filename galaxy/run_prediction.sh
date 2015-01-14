#!/bin/sh
# run_prediction $assembly $tag $strain $out $outfna
# Predicts where genes are in a genome assembly

ASSEMBLY=$1
TAG=$2
STRAIN=$3
OUT=$4
OUTFNA=$5

. `dirname $0`/cgpipelineGalaxyrc

if [ "" = "$OUT" ]; then
  echo "ERROR: need output parameter";
  echo "$0 ASSEMBLY.fasta TAG STRAIN OUT.gbk"
  exit 1;
fi;

TAG="$TAG"_

run_prediction $ASSEMBLY --tag_prefix=$TAG -strain_name=$STRAIN -crispr -o $OUT 2>&1
if [ $? -gt 0 ]; then
  echo "ERROR with prediction. See the log for details." 1>&2
  exit 1
fi
echo "Finished with gene prediction."
echo "Extracting all genes to nucleotide multi-fasta file..."
mv -v "$OUT.fna" "$OUTFNA"
if [ $? -gt 0 ]; then
  echo "ERROR with creating the fasta file. See the log for details." 1>&2
  exit 1
fi
