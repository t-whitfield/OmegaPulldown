#!/usr/bin/env python3
"""
Generate all bait x candidate pairs from FASTA files.

Reads two multi-sequence FASTA files (baits and candidates), extracts
protein IDs and sequences, and emits a TSV of all pairwise combinations.

The protein ID is taken from the FASTA header (first whitespace-delimited
token after '>').

Output TSV columns:
  pair_id, bait_id, bait_seq, candidate_id, candidate_seq

Pair IDs use the format: {bait}_with_{candidate}  (e.g. actin_with_HuR)
"""

import argparse
import sys
from itertools import product
from pathlib import Path


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


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--baits", required=True,
                        help="FASTA file containing bait sequences")
    parser.add_argument("--candidates", required=True,
                        help="FASTA file containing candidate sequences")
    parser.add_argument("--output", default="pairs.tsv", help="Output TSV file")
    parser.add_argument("--self-pairs", action="store_true",
                        help="Include self-pairs when a protein appears in both sets")
    args = parser.parse_args()

    baits = load_sequences_from_fasta(args.baits)
    candidates = load_sequences_from_fasta(args.candidates)

    if not baits:
        sys.exit(f"ERROR: no FASTA sequences found in {args.baits}")
    if not candidates:
        sys.exit(f"ERROR: no FASTA sequences found in {args.candidates}")

    print(f"Loaded {len(baits)} baits and {len(candidates)} candidates",
          file=sys.stderr)
    print(f"Generating {len(baits) * len(candidates)} pairs", file=sys.stderr)

    with open(args.output, "w") as out:
        out.write("pair_id\tbait_id\tbait_seq\tcandidate_id\tcandidate_seq\n")
        n = 0
        for bait_id, cand_id in product(sorted(baits.keys()),
                                         sorted(candidates.keys())):
            if not args.self_pairs and bait_id == cand_id:
                continue
            pair_id = f"{bait_id}_with_{cand_id}"
            out.write(f"{pair_id}\t{bait_id}\t{baits[bait_id]}"
                      f"\t{cand_id}\t{candidates[cand_id]}\n")
            n += 1

    print(f"Wrote {n} pairs to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
