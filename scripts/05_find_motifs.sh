#!/bin/bash
# Identify recognition motifs, twice, by different methods.
#
# MotifMaker searches the calls against the genome directly. MEME is
# run on sequence windows around the strongest calls as an independent
# check, and a third time on the 6mA calls that match none of the
# motifs already found, which is the test for an additional motif.

set -eu
source config.sh

conda activate smrt_env
motifMaker find -f "$WORK/assembly.fasta" -g "$WORK/basemods.gff" \
    -o "$WORK/motifs.csv"
cat "$WORK/motifs.csv"

mkdir -p "$WORK/meme"
conda activate ipd_env

# 41 bp windows centred on the strongest calls of each base type.
awk '!/^#/ && $3=="m6A"' "$WORK/basemods.gff" | sort -k6,6nr | head -1000 |
  grep -o 'context=[ACGT]*' | cut -d= -f2 |
  awk '{print ">site"NR"\n"$0}' > "$WORK/meme/top_6mA.fa"

awk '!/^#/ && $3=="m4C"' "$WORK/basemods.gff" | sort -k6,6nr | head -1000 |
  grep -o 'context=[ACGT]*' | cut -d= -f2 |
  awk '{print ">site"NR"\n"$0}' > "$WORK/meme/top_4mC.fa"

# Calls above the quality threshold matching none of the three 6mA
# motifs. The substring offsets pick the motif out of the 41 bp window.
awk '!/^#/ && $3=="m6A" && $6>=30' "$WORK/basemods.gff" |
  grep -o 'context=[ACGT]*' | cut -d= -f2 |
  awk '{cacag=substr($0,19,4); gtwwac=substr($0,17,6); typeI=substr($0,20,10);
        if (cacag!="ACAG" && gtwwac !~ /^GT[AT][AT]AC$/ && typeI !~ /^GAG.{3}CTC$/)
          print ">unassigned"NR"\n"$0}' > "$WORK/meme/unassigned_6mA.fa"

conda activate meme_env
meme "$WORK/meme/top_6mA.fa" -dna -oc "$WORK/meme/out_6mA" \
     -mod zoops -nmotifs 6 -minw 4 -maxw 14 -revcomp -markov_order 2
meme "$WORK/meme/top_4mC.fa" -dna -oc "$WORK/meme/out_4mC" \
     -mod zoops -nmotifs 4 -minw 4 -maxw 12 -revcomp -markov_order 2
meme "$WORK/meme/unassigned_6mA.fa" -dna -oc "$WORK/meme/out_unassigned" \
     -mod zoops -nmotifs 4 -minw 4 -maxw 14 -revcomp -markov_order 2

grep -hE "^MOTIF" "$WORK"/meme/out_*/meme.txt
