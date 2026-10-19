#!/bin/bash
# Annotate the assembly.
#
# Bakta is run on the assembly as produced by Flye. A second copy is
# reoriented to start at dnaA, which is the conventional presentation
# for a bacterial chromosome, but reorientation rotates and in one case
# reverse complements the sequence, so the two annotations share no
# coordinates. Everything in the paper uses the unreoriented frame;
# script 11 reconciles the two.

set -eu
source config.sh

conda activate bakta_env
bakta --db "$BAKTA_DB" --output "$WORK/bakta" --prefix wtBTH \
      --genus Burkholderia --species thailandensis --strain E264 \
      --threads "$THREADS" "$WORK/assembly.fasta"

conda activate dnaapler_env
dnaapler all --input "$WORK/assembly.fasta" --output "$WORK/dnaapler" \
             --prefix reoriented --threads "$THREADS"
