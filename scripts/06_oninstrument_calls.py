#!/usr/bin/env python3
"""Map the instrument's own modification calls onto the assembly.

The basecaller writes per-read modification probabilities into the MM
and ML tags of the unaligned BAM, using a model matched to the
sequencing chemistry. Those calls are a useful check on the kinetics
analysis, which used an older model, but they live in read coordinates
and have to be placed on the genome first.

Rather than realign, which would mean another format conversion and
risk losing the tags, this reads the tags from the unaligned file and
takes coordinates from the alignment produced in script 03, joining the
two on read name. Each call is required to land on its own canonical
base, which catches any coordinate error immediately.
"""

import argparse
import collections
import re
import sys

import pysam

CANONICAL = {"a": "A", "m": "C", "h": "C"}
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
MIN_PROBABILITY = 250  # out of 255, i.e. 0.98

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]",
         "S": "[GC]", "W": "[AT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]",
         "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "."}
IUPAC_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "R": "Y", "Y": "R",
                    "S": "S", "W": "W", "K": "M", "M": "K", "B": "V", "D": "H",
                    "H": "D", "V": "B", "N": "N"}

MOTIFS = [("CACAG", 3), ("GTWWAC", 4), ("GAGHNNNCTC", 1)]


def as_regex(motif):
    return "".join(IUPAC[base] for base in motif)


def reverse_complement(motif):
    return "".join(IUPAC_COMPLEMENT[base] for base in reversed(motif))


def read_fasta(path):
    sequences, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            sequences[name] = []
        else:
            sequences[name].append(line.strip().upper())
    return {name: "".join(parts) for name, parts in sequences.items()}


def motif_positions(sequences, motif, offset):
    """Genomic positions of the modified base, both strands."""
    forward, reverse = as_regex(motif), as_regex(reverse_complement(motif))
    positions = set()
    for name, sequence in sequences.items():
        for match in re.finditer("(?=(%s))" % forward, sequence):
            positions.add((name, match.start() + offset + 1))
        target = forward if reverse_complement(motif) == motif else reverse
        for match in re.finditer("(?=(%s))" % target, sequence):
            positions.add((name, match.start() + len(motif) - offset))
    return positions


def collect_read_calls(unaligned_bam):
    """Confident modification calls, keyed by read name."""
    calls = {}
    bam = pysam.AlignmentFile(unaligned_bam, "rb", check_sq=False)
    for read in bam:
        modifications = read.modified_bases
        if not modifications:
            continue
        confident = []
        for key, values in modifications.items():
            code = key[2] if isinstance(key, tuple) and len(key) > 2 else str(key)
            if code not in CANONICAL:
                continue
            confident += [(position, code) for position, probability in values
                          if probability >= MIN_PROBABILITY]
        if confident:
            calls[read.query_name] = (len(read.query_sequence), confident)
    return calls


def place_on_genome(aligned_bam, sequences, calls):
    """Translate read coordinates to genome coordinates."""
    placed = collections.defaultdict(collections.Counter)
    matched = mismatched = 0
    bam = pysam.AlignmentFile(aligned_bam, "rb")
    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        entry = calls.get(read.query_name)
        if not entry:
            continue
        length, confident = entry
        pairs = dict(read.get_aligned_pairs(matches_only=True))
        for query_position, code in confident:
            # The stored positions index the original read, so reverse
            # them for reads aligned to the minus strand.
            index = length - 1 - query_position if read.is_reverse else query_position
            reference_position = pairs.get(index)
            if reference_position is None:
                continue
            base = sequences[read.reference_name][reference_position]
            expected = CANONICAL[code]
            if base == expected:
                strand = "+"
            elif base == COMPLEMENT.get(expected):
                strand = "-"
            else:
                mismatched += 1
                continue
            matched += 1
            placed[code][(read.reference_name, reference_position + 1)] += 1
    return placed, matched, mismatched


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("assembly")
    parser.add_argument("unaligned_bam")
    parser.add_argument("aligned_bam")
    parser.add_argument("output")
    args = parser.parse_args()

    sequences = read_fasta(args.assembly)
    calls = collect_read_calls(args.unaligned_bam)
    print(f"reads carrying a confident call: {len(calls):,}")

    placed, matched, mismatched = place_on_genome(
        args.aligned_bam, sequences, calls)
    total = matched + mismatched
    print(f"calls placed on the genome: {matched:,} of {total:,} "
          f"({100 * matched / total:.1f}% on the expected base)")
    if total and matched / total < 0.9:
        sys.exit("too many calls fall on the wrong base; check the coordinate "
                 "convention before using these results")

    adenine_calls = set(placed.get("a", {}))
    with open(args.output, "w") as out:
        out.write("motif\tgenomic_positions\trecovered\tpercent\n")
        for motif, offset in MOTIFS:
            positions = motif_positions(sequences, motif, offset)
            recovered = len(positions & adenine_calls)
            out.write(f"{motif}\t{len(positions)}\t{recovered}\t"
                      f"{100 * recovered / len(positions):.1f}\n")
            print(f"  {motif:12s} {recovered:,} of {len(positions):,} "
                  f"({100 * recovered / len(positions):.1f}%)")

    cytosine_calls = placed.get("m", {})
    print(f"\n5-methylcytosine calls: {len(cytosine_calls):,} positions")
    print("These are CpG context only and are assessed in the paper, "
          "Technical Validation.")


if __name__ == "__main__":
    main()
