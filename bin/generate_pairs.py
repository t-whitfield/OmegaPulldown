#!/usr/bin/env python3
"""
Generate the pair list that drives the rest of the pipeline.

Two input modes, producing an identical output TSV so nothing downstream
changes:

  1. CROSS-PRODUCT (--baits/--candidates): every bait x candidate pair.
     The classic pulldown design -- one or a few baits against a candidate
     list, e.g. the hits from one AP-MS experiment.

  2. PAIR LIST (--pairs-tsv/--sequences): only the pairs named in a TSV.
     For adding structure-based ranking to an existing edge list, e.g.
     interacting pairs reported by proteomics.

The protein ID is the first whitespace-delimited token of the FASTA header.

Output TSV columns:
  pair_id, bait_id, bait_seq, candidate_id, candidate_seq

Pair IDs use the format: {bait}_with_{candidate}  (e.g. actin_with_HuR)

UNORDERED PAIRS
---------------
A complex is unordered: A_with_B and B_with_A are the same system, and
predicting both wastes half the GPU budget. Both modes therefore collapse
pairs on the canonical key tuple(sorted((a, b))) and keep the first
orientation seen. This is a no-op when baits and candidates are disjoint
(the usual pulldown case) and halves the work exactly when they overlap --
e.g. an all-vs-all run with --baits and --candidates set to the same file,
which previously emitted N^2 - N ordered pairs instead of N(N-1)/2.

Which protein ends up as the "bait" of a collapsed pair is determined by
iteration order (sorted, so deterministic and reproducible). That choice is
purely cosmetic: it fixes chain order (bait is chain A) and the pair_id,
nothing else.
"""

import argparse
import csv
import sys
from itertools import product


def parse_fasta(fasta_path):
    """Yield (id, sequence) tuples from a FASTA file."""
    header = None
    seq_lines = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line[1:].split()[0]
                seq_lines = []
            elif line:
                seq_lines.append(line)
    if header is not None:
        yield header, "".join(seq_lines)


def load_sequences_from_fasta(fasta_path):
    """Load all sequences from a single multi-sequence FASTA file."""
    seqs = {}
    for prot_id, seq in parse_fasta(fasta_path):
        if prot_id in seqs:
            print(f"WARNING: duplicate ID '{prot_id}' in {fasta_path}, skipping",
                  file=sys.stderr)
            continue
        seqs[prot_id] = seq
    return seqs


def read_pair_list(pairs_path, seqs):
    """
    Read (id_a, id_b) rows from a TSV or CSV. Only the first two columns are
    used; any further columns (MS score, FDR, spectral counts) are ignored
    here: they are not threaded through to the rankings yet.

    A header row is detected rather than assumed: if neither of the first
    row's two IDs resolves against the FASTA, it is treated as a header.
    That way both a bare two-column list and one exported with column names
    work without a flag.
    """
    rows = []
    with open(pairs_path, newline="") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        delim = "\t" if "\t" in sample.splitlines()[0] else ","
        for lineno, row in enumerate(csv.reader(fh, delimiter=delim), start=1):
            row = [c.strip() for c in row if c.strip() != ""]
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 2:
                sys.exit(f"ERROR: {pairs_path} line {lineno}: need two protein "
                         f"IDs, got {row!r}")
            rows.append((row[0], row[1], lineno))

    if not rows:
        sys.exit(f"ERROR: no pairs found in {pairs_path}")

    a, b, _ = rows[0]
    if a not in seqs and b not in seqs:
        print(f"Treating first line of {pairs_path} as a header: {a!r}, {b!r}",
              file=sys.stderr)
        rows = rows[1:]
        if not rows:
            sys.exit(f"ERROR: {pairs_path} contains only a header")
    return rows


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--baits", help="FASTA of bait sequences (cross-product mode)")
    parser.add_argument("--candidates", help="FASTA of candidate sequences (cross-product mode)")
    parser.add_argument("--pairs-tsv", dest="pairs_tsv",
                        help="TSV/CSV whose first two columns name the pairs to predict (pair-list mode)")
    parser.add_argument("--sequences",
                        help="FASTA holding every protein referenced by --pairs-tsv")
    parser.add_argument("--output", default="pairs.tsv", help="Output TSV file")
    parser.add_argument("--self-pairs", action="store_true",
                        help="Cross-product mode only: include homodimers (A with A). "
                             "In pair-list mode an explicitly listed self-pair is always "
                             "kept, since the list states the user's intent.")
    args = parser.parse_args()

    pair_mode = bool(args.pairs_tsv)
    if pair_mode:
        if not args.sequences:
            sys.exit("ERROR: --pairs-tsv requires --sequences (a FASTA with every protein it names)")
        if args.baits or args.candidates:
            sys.exit("ERROR: --pairs-tsv is mutually exclusive with --baits/--candidates")
    else:
        if not (args.baits and args.candidates):
            sys.exit("ERROR: need either --baits and --candidates, or --pairs-tsv and --sequences")

    seen = set()
    dupes = 0
    n = 0

    with open(args.output, "w") as out:
        out.write("pair_id\tbait_id\tbait_seq\tcandidate_id\tcandidate_seq\n")

        if pair_mode:
            seqs = load_sequences_from_fasta(args.sequences)
            if not seqs:
                sys.exit(f"ERROR: no FASTA sequences found in {args.sequences}")
            rows = read_pair_list(args.pairs_tsv, seqs)
            print(f"Loaded {len(seqs)} sequences and {len(rows)} listed pairs",
                  file=sys.stderr)

            missing = {}
            for a, b, lineno in rows:
                for pid in (a, b):
                    if pid not in seqs:
                        missing.setdefault(pid, lineno)
            if missing:
                shown = list(missing.items())[:10]
                detail = "\n".join(f"    {pid}  (first seen line {ln})"
                                   for pid, ln in shown)
                more = "" if len(missing) <= 10 else f"\n    ... and {len(missing) - 10} more"
                sys.exit(f"ERROR: {len(missing)} protein ID(s) in {args.pairs_tsv} "
                         f"are absent from {args.sequences}:\n{detail}{more}")

            for a, b, _ in rows:
                key = tuple(sorted((a, b)))
                if key in seen:
                    dupes += 1
                    continue
                seen.add(key)
                out.write(f"{a}_with_{b}\t{a}\t{seqs[a]}\t{b}\t{seqs[b]}\n")
                n += 1
        else:
            baits = load_sequences_from_fasta(args.baits)
            candidates = load_sequences_from_fasta(args.candidates)
            if not baits:
                sys.exit(f"ERROR: no FASTA sequences found in {args.baits}")
            if not candidates:
                sys.exit(f"ERROR: no FASTA sequences found in {args.candidates}")
            print(f"Loaded {len(baits)} baits and {len(candidates)} candidates",
                  file=sys.stderr)

            for bait_id, cand_id in product(sorted(baits.keys()),
                                            sorted(candidates.keys())):
                if not args.self_pairs and bait_id == cand_id:
                    continue
                key = tuple(sorted((bait_id, cand_id)))
                if key in seen:
                    dupes += 1
                    continue
                seen.add(key)
                out.write(f"{bait_id}_with_{cand_id}\t{bait_id}\t{baits[bait_id]}"
                          f"\t{cand_id}\t{candidates[cand_id]}\n")
                n += 1

    if dupes:
        print(f"Collapsed {dupes} duplicate unordered pair(s) "
              f"(A_with_B == B_with_A)", file=sys.stderr)
    print(f"Wrote {n} pairs to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
