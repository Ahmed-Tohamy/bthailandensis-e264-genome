#!/bin/bash
# Call modified bases from polymerase kinetics.
#
# Kinetics are measured per strand, but a HiFi read is a consensus of
# both, carrying forward and reverse measurements in separate tags on
# one record. bystrandify splits each read into two pseudo-reads so
# that kineticsTools sees the per-strand representation it expects.
# The header is restored afterwards because bystrandify rewrites the
# read group, which pbmm2 then rejects.
#
# The kineticsTools release used here predates the Revio platform and
# ships no model for P2-C2 chemistry. With automatic model selection
# it exits reporting that the chemistry cannot be identified, so a
# Sequel II model is given explicitly. Kinetic values produced this way
# are not calibrated in absolute terms; see the paper, Technical
# Validation, for the two checks applied in consequence.

set -eu
source config.sh

conda activate smrt_env
samtools view -H "$WORK/merged_unaligned.bam" > "$WORK/original_header.sam"

"$SMRTLINK/ccs-kinetics-bystrandify" "$WORK/merged_unaligned.bam" \
    "$WORK/bystrand.bam"
samtools reheader "$WORK/original_header.sam" "$WORK/bystrand.bam" \
    > "$WORK/bystrand_fixed.bam"

pbmm2 align "$WORK/assembly.fasta" "$WORK/bystrand_fixed.bam" \
    "$WORK/bystrand_aligned.bam" --sort -j "$THREADS" --preset HIFI
pbindex "$WORK/bystrand_aligned.bam"

ipdSummary "$WORK/bystrand_aligned.bam" \
    --reference "$WORK/assembly.fasta" \
    --gff "$WORK/basemods.gff" --csv "$WORK/basemods.csv" \
    --ipdModel "$KINETICS_MODEL" \
    --identify m6A,m4C --methylFraction \
    --numWorkers "$THREADS" --minCoverage 5

awk '!/^#/{print $3}' "$WORK/basemods.gff" | sort | uniq -c
