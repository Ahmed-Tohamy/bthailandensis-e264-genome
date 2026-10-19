#!/bin/bash
# Align the assembly against the three published assemblies of E264.
#
# Whether a whole replicon is stored forward or reverse complemented is
# an assembly convention, not a property of the organism, so a query
# replicon is flipped when most of its aligned length comes back
# inverted. Without that step every comparison looks rearranged.

set -eu
source config.sh

mkdir -p "$WORK/compare"
conda activate amplicon_env
for ACCESSION in GCF_000012365.1 GCA_003568605.1 GCA_000152285.1; do
  if [ ! -s "$WORK/compare/$ACCESSION.fna" ]; then
    datasets download genome accession "$ACCESSION" \
      --filename "$WORK/compare/$ACCESSION.zip"
    unzip -o -j "$WORK/compare/$ACCESSION.zip" "*genomic.fna" \
      -d "$WORK/compare" && mv "$WORK"/compare/*genomic.fna \
      "$WORK/compare/$ACCESSION.fna"
  fi
done

conda activate mummer_env
for ACCESSION in GCF_000012365.1 GCA_003568605.1 GCA_000152285.1; do
  REFERENCE="$WORK/compare/$ACCESSION.fna"
  nucmer --maxmatch -c 2000 -p "$WORK/compare/synteny_$ACCESSION" \
         "$REFERENCE" "$WORK/assembly.fasta"
  delta-filter -1 -l 20000 "$WORK/compare/synteny_$ACCESSION.delta" \
    > "$WORK/compare/synteny_$ACCESSION.filtered.delta"
  show-coords -rclTH "$WORK/compare/synteny_$ACCESSION.filtered.delta" \
    > "$WORK/compare/synteny_$ACCESSION.coords"

  # Finer anchoring for single base differences.
  nucmer --maxmatch -c 200 -p "$WORK/compare/variants_$ACCESSION" \
         "$REFERENCE" "$WORK/assembly.fasta"
  delta-filter -1 -l 1000 "$WORK/compare/variants_$ACCESSION.delta" \
    > "$WORK/compare/variants_$ACCESSION.filtered.delta"
  show-snps -ClrTH "$WORK/compare/variants_$ACCESSION.filtered.delta" \
    > "$WORK/compare/variants_$ACCESSION.snps"

  echo "$ACCESSION: $(wc -l < "$WORK/compare/synteny_$ACCESSION.coords") blocks, \
$(wc -l < "$WORK/compare/variants_$ACCESSION.snps") differences"
done

conda activate fastani_env
ls "$WORK"/compare/GC*.fna > "$WORK/compare/reference_list.txt"
fastANI -q "$WORK/assembly.fasta" --rl "$WORK/compare/reference_list.txt" \
        -o "$WORK/compare/ani.tsv" -t "$THREADS"
cat "$WORK/compare/ani.tsv"
