#!/usr/bin/env python3
"""Decide which base differences are real and which are reference errors.

A difference between one genome and one reference cannot distinguish a
true variant from an error in that reference. Calling against two
assemblies of the same strain separates them: a difference seen against
both is a property of this isolate, one seen against only one is more
likely a problem in that assembly.

Differences cannot be matched on reference coordinates, because the
alignment blocks differ between comparisons and the same difference can
fall in different blocks. They are matched instead on the sequence
surrounding them in our own assembly, which is the same object in every
comparison.
"""

import argparse
import collections
import os

CONTEXT = 10  # bases either side


def read_fasta(path):
    sequences, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            sequences[name] = []
        else:
            sequences[name].append(line.strip().upper())
    return {name: "".join(parts) for name, parts in sequences.items()}


def read_differences(path, sequences):
    """Key each difference by its context in our assembly."""
    differences = {}
    if not os.path.exists(path):
        return differences
    for line in open(path):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 5:
            continue
        contig = next((f for f in fields if f in sequences), None)
        if contig is None:
            continue
        numbers = [f for f in fields[:6] if f.lstrip("-").isdigit()]
        if len(numbers) < 2:
            continue
        position = int(numbers[1])
        sequence = sequences[contig]
        if position <= CONTEXT or position + CONTEXT > len(sequence):
            continue
        context = sequence[position - CONTEXT - 1:position + CONTEXT]
        differences[(contig, context)] = (fields[1], fields[2], position)
    return differences


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("assembly")
    parser.add_argument("--snps", nargs="+", required=True,
                        metavar="LABEL=PATH",
                        help="one show-snps output per reference")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    sequences = read_fasta(args.assembly)
    per_reference = {}
    for item in args.snps:
        label, path = item.split("=", 1)
        per_reference[label] = read_differences(path, sequences)
        print(f"{label}: {len(per_reference[label]):,} differences")

    every_key = set().union(*per_reference.values())
    labels = list(per_reference)

    rows = []
    for key in every_key:
        seen_in = [l for l in labels if key in per_reference[l]]
        reference_base, our_base, position = per_reference[seen_in[0]][key]
        rows.append((key[0], position, reference_base, our_base, seen_in))
    rows.sort(key=lambda row: (row[0], row[1]))

    tally = collections.Counter(tuple(row[4]) for row in rows)
    print("\ndifferences by the references supporting them")
    for combination, count in sorted(tally.items(), key=lambda x: -x[1]):
        print(f"  {', '.join(combination):40s} {count:,}")

    with open(args.output, "w") as out:
        out.write("replicon\tposition\treference_base\tour_base\t"
                  + "\t".join(f"differs_from_{l}" for l in labels)
                  + "\tn_references\n")
        for replicon, position, reference_base, our_base, seen_in in rows:
            flags = "\t".join("yes" if l in seen_in else "no" for l in labels)
            out.write(f"{replicon}\t{position}\t{reference_base}\t{our_base}\t"
                      f"{flags}\t{len(seen_in)}\n")
    print(f"\nwritten to {args.output}")


if __name__ == "__main__":
    main()
