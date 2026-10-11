#!/usr/bin/env python3
"""
Convert an a3m MSA file to Boltz-2 CSV format.

Boltz-2 expects MSAs as CSV files with two columns:
  - sequence: the aligned protein sequence
  - key:      a unique identifier used to pair rows across chains

For paired MSAs, sequences from the same organism (same TaxID) share
the same key across different chain CSV files. This enables Boltz-2 to
perform evolutionary covariance analysis across interacting chains.

The script extracts TaxID from UniRef-style headers:
  >UniRef100_A0A0A0MQR9 description n=1 Tax=Homo sapiens TaxID=9606 ...

Sequences without a parseable TaxID receive key -1 (unpaired).
The query sequence (first entry) always receives key 0.

Input:  a3m file (from mmseqs2, colabfold_search, or jackhmmer)
Output: CSV file suitable for Boltz-2 YAML msa field
"""

import argparse
import csv
import re
import sys


def parse_a3m(a3m_path):
    """
    Parse an a3m file, yielding (header, sequence) tuples.
    a3m format: lowercase letters are insertions (removed for alignment).
    """
    header = None
    seq_lines = []
    with open(a3m_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    raw_seq = "".join(seq_lines)
                    # Remove lowercase insertion characters (a3m convention)
                    aligned_seq = re.sub(r"[a-z]", "", raw_seq)
                    yield header, aligned_seq
                header = line[1:]
                seq_lines = []
            elif line.startswith("#"):
                continue  # skip comment lines
            elif line:
                seq_lines.append(line)
    if header is not None:
        raw_seq = "".join(seq_lines)
        aligned_seq = re.sub(r"[a-z]", "", raw_seq)
        yield header, aligned_seq


def extract_taxid(header):
    """
    Extract TaxID from a UniRef-style header.
    Returns the TaxID string, or None if not found.
    """
    # Pattern: TaxID=12345 (possibly followed by space or end of string)
    match = re.search(r"TaxID=(\d+)", header)
    if match:
        return match.group(1)

    # Also try OX= format (UniProt style)
    match = re.search(r"OX=(\d+)", header)
    if match:
        return match.group(1)

    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("a3m", help="Input a3m file")
    parser.add_argument("--output", "-o", required=True, help="Output CSV file")
    parser.add_argument("--max-seqs", type=int, default=8192,
                        help="Maximum number of MSA sequences to retain (default: 8192)")
    args = parser.parse_args()

    rows_written = 0

    with open(args.output, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["sequence", "key"])

        for i, (header, seq) in enumerate(parse_a3m(args.a3m)):
            if rows_written >= args.max_seqs:
                break

            if i == 0:
                # Query sequence always gets key 0
                key = 0
            else:
                taxid = extract_taxid(header)
                if taxid is not None:
                    key = int(taxid)
                else:
                    # Unpaired sequences get key -1
                    key = -1

            writer.writerow([seq, key])
            rows_written += 1

    print(f"Converted {rows_written} sequences from {args.a3m} -> {args.output}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
