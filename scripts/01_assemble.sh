#!/bin/bash
# Assemble the HiFi reads into two closed chromosomes.
#
# The two SMRT Cells are merged as unaligned BAM rather than through
# FASTQ, because FASTQ conversion discards the per-base kinetic tags
# that the modification calling in script 04 depends on. The FASTQ
# written here is used only by Flye, which ignores the tags anyway.

set -eu
source config.sh

conda activate smrt_env
samtools merge -@ "$THREADS" -f "$WORK/merged_unaligned.bam" \
  "$READS_CELL1" "$READS_CELL2"
samtools fastq -@ "$THREADS" "$WORK/merged_unaligned.bam" > "$WORK/reads.fq"

conda activate flye_env
flye --pacbio-hifi "$WORK/reads.fq" --out-dir "$WORK/flye" --threads "$THREADS"

# Flye names contigs in the order it finishes them, so the larger
# replicon is contig_2 here. Rename to the convention used in the
# paper, where chromosome 1 is the larger of the two.
conda activate amplicon_env
seqkit replace -p '^contig_2$' -r 'chromosome_1' "$WORK/flye/assembly.fasta" |
  seqkit replace -p '^contig_1$' -r 'chromosome_2' > "$WORK/assembly.fasta"
seqkit stats -a "$WORK/assembly.fasta"
