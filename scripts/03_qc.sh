#!/bin/bash
# Assess assembly completeness and map the reads back.
#
# The depth file written at the end is used by scripts 05 and 10, so
# it covers every position including those with no coverage.

set -eu
source config.sh

conda activate checkm2_env
checkm2 predict --input "$WORK/assembly.fasta" --output-directory "$WORK/checkm2" \
                --threads "$THREADS" -x fasta

conda activate busco5_env
busco --in "$WORK/assembly.fasta" --out busco --out_path "$WORK" \
      --mode genome --lineage_dataset burkholderiales_odb10 --cpu "$THREADS"

conda activate assembly_preprocess_env
minimap2 -ax map-hifi -t "$THREADS" "$WORK/assembly.fasta" "$WORK/reads.fq" \
  > "$WORK/aligned.sam"

conda activate smrt_env
samtools sort -@ "$THREADS" -o "$WORK/aligned.bam" "$WORK/aligned.sam"
samtools index "$WORK/aligned.bam"
rm "$WORK/aligned.sam"
samtools flagstat "$WORK/aligned.bam"
samtools depth -a "$WORK/aligned.bam" > "$WORK/depth.txt"
