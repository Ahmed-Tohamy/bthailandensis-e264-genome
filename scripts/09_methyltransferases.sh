#!/bin/bash
# Identify methyltransferases and ask how widely they are conserved.
#
# REBASE is the reference collection for restriction and modification
# proteins. Hits give a candidate enzyme for each motif, though an
# assignment made this way is a prediction: confirming it would mean
# deleting the enzyme and resequencing.
#
# The conservation step asks whether an enzyme is present across the
# species or confined to this strain, which tells a reader whether to
# expect the same motifs elsewhere.

set -eu
source config.sh

mkdir -p "$WORK/rebase"
if [ ! -s "$WORK/rebase/proteins.faa" ]; then
  curl -s "ftp://ftp.neb.com/pub/rebase/protein_org_seqs.txt" \
    -o "$WORK/rebase/raw.txt"
  # The distributed file carries a text preamble and wrapped lines that
  # makeblastdb rejects, so keep only well formed records.
  awk '/^>/{header=$0; sequence=""; next}
       header && /^[ACDEFGHIKLMNPQRSTVWYXBZJUO*]+$/{sequence=sequence $0; next}
       header && sequence {print header"\n"sequence; header=""; sequence=""}' \
    "$WORK/rebase/raw.txt" > "$WORK/rebase/proteins.faa"
fi
echo "REBASE proteins: $(grep -c '^>' "$WORK/rebase/proteins.faa")"

conda activate abricate_env
makeblastdb -in "$WORK/rebase/proteins.faa" -dbtype prot \
            -out "$WORK/rebase/db" > /dev/null
blastp -query "$WORK/bakta/wtBTH.faa" -db "$WORK/rebase/db" \
       -outfmt '6 qseqid sseqid pident length qlen evalue bitscore stitle' \
       -max_target_seqs 3 -evalue 1e-30 -qcov_hsp_perc 50 \
       -num_threads "$THREADS" > "$WORK/rebase/hits.tsv"
echo "proteins with a REBASE match: $(cut -f1 "$WORK/rebase/hits.tsv" | sort -u | wc -l)"

# Conservation across the species. Accessions are those compared in the
# paper; edit the list to compare against others.
mkdir -p "$WORK/conservation"
conda activate amplicon_env
# E444, FDAARGOS_236, FDAARGOS_425, AW34-19p, 2002721643 and BPM.
# Two of these, FDAARGOS_236 and FDAARGOS_425, turn out to be further
# deposits of E264 rather than distinct strains; the ANI reported by
# script 07 is what shows this, and the paper reads them that way.
for ACCESSION in GCA_000567945.1 GCA_002891115.1 GCA_002588555.1 \
                 GCA_016774495.1 GCA_001718635.1 GCA_011578485.1; do
  if [ ! -s "$WORK/conservation/$ACCESSION.faa" ]; then
    datasets download genome accession "$ACCESSION" --include protein \
      --filename "$WORK/conservation/$ACCESSION.zip" || continue
    unzip -o -j "$WORK/conservation/$ACCESSION.zip" "*protein.faa" \
      -d "$WORK/conservation" && \
      mv "$WORK"/conservation/protein.faa "$WORK/conservation/$ACCESSION.faa"
  fi
done

conda activate abricate_env
for FAA in "$WORK"/conservation/*.faa; do
  LABEL=$(basename "$FAA" .faa)
  makeblastdb -in "$FAA" -dbtype prot -out "$WORK/conservation/db_$LABEL" > /dev/null
  blastp -query "$WORK/methyltransferases.faa" -db "$WORK/conservation/db_$LABEL" \
         -outfmt '6 qseqid pident length qlen bitscore' -max_target_seqs 1 \
         -evalue 1e-20 -num_threads "$THREADS" |
    sort -k1,1 -k5,5gr | awk -v s="$LABEL" '!seen[$1]++{print $1"\t"s"\t"$2}'
done > "$WORK/conservation/identity.tsv"
cat "$WORK/conservation/identity.tsv"
