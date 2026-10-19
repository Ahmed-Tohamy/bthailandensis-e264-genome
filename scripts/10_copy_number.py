#!/usr/bin/env python3
"""Measure the copy number of an amplified segment from read depth.

Depth is a usable proxy for copy number only after two corrections.
HiFi coverage varies with local GC content, so window depths are
compared against a depth-versus-GC curve built from the rest of the
genome rather than against a flat average. And depth varies from window
to window simply because reads are sampled, so control regions of the
same size elsewhere in the genome give the spread expected when copy
number is uniform.

Tandem arrangement is established separately, from reads whose two
alignments run from the end of the segment back to its start.
"""

import argparse
import collections
import re
import statistics

import pysam

WINDOW = 10000


def read_fasta(path):
    sequences, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            sequences[name] = []
        else:
            sequences[name].append(line.strip().upper())
    return {name: "".join(parts) for name, parts in sequences.items()}


def window_depths(depth_file, replicon):
    totals = collections.Counter()
    counts = collections.Counter()
    for line in open(depth_file):
        name, position, depth = line.split()
        if name != replicon:
            continue
        window = int(position) // WINDOW
        totals[window] += int(depth)
        counts[window] += 1
    return {w: totals[w] / counts[w] for w in totals if counts[w] > WINDOW * 0.9}


def gc_fraction(sequence):
    return 100 * (sequence.count("G") + sequence.count("C")) / len(sequence)


def reference_length(cigar):
    return sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
               if op in "MDN=X")


def junction_at_position(bam_path, replicon, start, end, position):
    """Reads spanning a position, and how many take the tandem junction.

    Counted at a position rather than across the whole segment, because
    the informative quantity is the fraction of molecules at the end of
    the segment that continue into another copy rather than into
    downstream sequence. That fraction is what says the culture holds
    both arrangements.
    """
    spanning = junction = 0
    bam = pysam.AlignmentFile(bam_path, "rb")
    for read in bam.fetch(replicon, position - 1, position):
        if read.is_unmapped or read.is_secondary:
            continue
        spanning += 1
        if not read.has_tag("SA"):
            continue
        for entry in read.get_tag("SA").split(";"):
            if not entry:
                continue
            fields = entry.split(",")
            if len(fields) < 4 or fields[0] != replicon:
                continue
            other_start = int(fields[1])
            # A supplementary alignment landing back near the start of
            # the segment means this molecule carries the next copy.
            if abs(other_start - start) < 20000:
                junction += 1
                break
    return spanning, junction


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("assembly")
    parser.add_argument("depth_file")
    parser.add_argument("aligned_bam")
    parser.add_argument("--replicon", required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    sequences = read_fasta(args.assembly)
    sequence = sequences[args.replicon]
    depths = window_depths(args.depth_file, args.replicon)

    inside, outside = {}, {}
    for window, depth in depths.items():
        start = window * WINDOW
        target = inside if args.start <= start and start + WINDOW <= args.end else outside
        target[window] = depth

    # Expected depth as a function of GC, from outside the segment.
    by_gc = collections.defaultdict(list)
    for window, depth in outside.items():
        by_gc[round(gc_fraction(sequence[window * WINDOW:(window + 1) * WINDOW]))].append(depth)
    curve = {gc: statistics.median(values)
             for gc, values in by_gc.items() if len(values) >= 5}

    def expected(gc):
        gc = round(gc)
        for offset in range(10):
            for candidate in (gc - offset, gc + offset):
                if candidate in curve:
                    return curve[candidate]
        return statistics.median(curve.values())

    # Compare against the flanking sequence rather than the whole
    # replicon. An amplification is a local event, and depth drifts
    # across a chromosome for reasons that have nothing to do with it.
    flank = (args.end - args.start)
    flanking = {w: d for w, d in outside.items()
                if args.start - flank <= w * WINDOW <= args.end + flank}
    background = statistics.mean(flanking.values()) if flanking \
        else statistics.mean(outside.values())
    print(f"flanking windows used as background: {len(flanking)}")
    raw, corrected = [], []
    rows = []
    for window in sorted(inside):
        start = window * WINDOW
        gc = gc_fraction(sequence[start:start + WINDOW])
        raw.append(inside[window] / background)
        corrected.append(inside[window] / expected(gc))
        rows.append((start, gc, inside[window], raw[-1], corrected[-1]))

    print(f"windows inside the segment: {len(inside)}")
    print(f"background depth: {background:.1f}")
    # Report the ratio of mean depths, not the median of per-window
    # ratios. The two differ when depth varies within the segment, and
    # the first is the quantity meant by copy number.
    inside_mean = statistics.mean(inside.values())
    print(f"raw ratio:       {inside_mean / background:.2f}  "
          f"per window {min(raw):.2f} to {max(raw):.2f}")
    print(f"GC corrected:    {statistics.median(corrected):.2f}  "
          f"per window {min(corrected):.2f} to {max(corrected):.2f}")

    print("\nreads at the end of the segment, and how many continue")
    print("into another copy rather than into downstream sequence:")
    # Sampled just past the end of the segment, where a molecule has
    # already committed to one arrangement or the other. Inside the
    # segment the question has not yet been put to the read.
    for offset in (-1000, 0, 1000, 2000):
        position = args.end + offset
        spanning, junction = junction_at_position(
            args.aligned_bam, args.replicon, args.start, args.end, position)
        print(f"  {args.replicon}:{position:,}  reads {spanning:>3}  "
              f"taking the junction {junction:>3}")

    with open(args.output, "w") as out:
        out.write("window_start\tgc_percent\tdepth\traw_ratio\tcorrected_ratio\n")
        for start, gc, depth, raw_ratio, corrected_ratio in rows:
            out.write(f"{start}\t{gc:.1f}\t{depth:.1f}\t{raw_ratio:.2f}\t"
                      f"{corrected_ratio:.2f}\n")
    print(f"written to {args.output}")


if __name__ == "__main__":
    main()
