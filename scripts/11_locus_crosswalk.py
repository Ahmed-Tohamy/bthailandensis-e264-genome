#!/usr/bin/env python3
"""Match our locus tags to those of the 2005 E264 reference.

Three sets of locus tags are in play: ours, the RefSeq tags of the
reference, and the original BTH tags that the literature uses. The
original tags survive in the reference annotation as old_locus_tag
attributes on gene features, while protein identifiers appear on the
CDS features beneath them, so the two have to be joined through the
parent relationship rather than read from one line.

Matching is done on protein sequence, which is unaffected by the
rotation and reverse complementation that distinguish our assembly's
coordinate frame from the reference.
"""

import argparse
import collections
import csv
import re


def attribute(text, key):
    match = re.search(rf"(?:^|;){key}=([^;]+)", text)
    return match.group(1) if match else None


def read_reference_annotation(gff_path):
    """protein accession to (RefSeq tag, original BTH tag, replicon)."""
    gene_to_original, tag_to_original = {}, {}
    protein_to_gene, protein_to_tag, protein_to_replicon = {}, {}, {}

    for line in open(gff_path):
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        attributes = fields[8]
        if fields[2] == "gene":
            identifier = attribute(attributes, "ID")
            original = attribute(attributes, "old_locus_tag")
            tag = attribute(attributes, "locus_tag")
            if identifier and original:
                gene_to_original[identifier] = original
            if tag and original:
                tag_to_original[tag] = original
        elif fields[2] == "CDS":
            protein = attribute(attributes, "protein_id")
            if not protein:
                continue
            protein_to_gene[protein] = attribute(attributes, "Parent")
            protein_to_tag[protein] = attribute(attributes, "locus_tag")
            protein_to_replicon[protein] = fields[0]

    lookup = {}
    for protein, gene in protein_to_gene.items():
        refseq_tag = protein_to_tag.get(protein, "")
        original = (gene_to_original.get(gene)
                    or tag_to_original.get(refseq_tag)
                    or "")
        lookup[protein] = (refseq_tag, original, protein_to_replicon.get(protein, ""))
    return lookup


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("blast_hits", help="blastp against the reference proteome")
    parser.add_argument("reference_gff")
    parser.add_argument("output")
    args = parser.parse_args()

    lookup = read_reference_annotation(args.reference_gff)
    print(f"reference proteins with annotation: {len(lookup):,}")
    with_original = sum(1 for v in lookup.values() if v[1])
    print(f"of which carry an original BTH tag: {with_original:,}")

    best = {}
    for row in csv.reader(open(args.blast_hits), delimiter="\t"):
        if len(row) < 7:
            continue
        query, score = row[0], float(row[6])
        if query not in best or score > best[query][0]:
            best[query] = (score, row)

    shared = collections.Counter(best[q][1][1] for q in best)

    with open(args.output, "w") as out:
        out.write("our_locus\tE264_locus_tag\tRefSeq_locus_tag\t"
                  "protein_accession\treplicon\tpercent_identity\tnote\n")
        for query in sorted(best):
            row = best[query][1]
            protein = row[1]
            refseq_tag, original, replicon = lookup.get(protein, ("", "", ""))
            identity = float(row[2])
            notes = []
            if identity < 99:
                notes.append("partial match")
            if shared[protein] > 1:
                notes.append("shares a reference protein with another query")
            out.write(f"{query}\t{original or 'not assigned'}\t{refseq_tag}\t"
                      f"{protein}\t{replicon}\t{identity:.1f}\t"
                      f"{'; '.join(notes)}\n")
    print(f"written to {args.output}")


if __name__ == "__main__":
    main()
