#!/bin/bash
# Pull out the methyltransferases whose motifs the paper reports.
#
# Script 09 compares these four against other genomes of the species.
# They are identified by the REBASE search in the first half of 09;
# the locus tags below are those the search returned for this assembly
# and will differ in any reannotation, so check them against
# tables/motif_assignments.tsv before rerunning.

set -eu
source config.sh

cat > "$WORK/methyltransferase_loci.txt" << 'LOCI'
FAADHO_05464
FAADHO_01359
FAADHO_04242
FAADHO_01013
LOCI

conda activate amplicon_env
seqkit grep -f "$WORK/methyltransferase_loci.txt" "$WORK/bakta/wtBTH.faa" \
  > "$WORK/methyltransferases.faa"

echo "selected $(grep -c '^>' "$WORK/methyltransferases.faa") of 4 proteins"
grep '^>' "$WORK/methyltransferases.faa" | cut -c1-70
