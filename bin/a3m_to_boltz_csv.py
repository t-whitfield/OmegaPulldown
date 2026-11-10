#!/usr/bin/env python3
"""
Convert an a3m MSA file to Boltz-2 CSV format.

Boltz-2 expects MSAs as CSV files with two columns:
  - sequence: the aligned protein sequence
  - key:      a unique identifier used to pair rows across chains

For paired MSAs, sequences from the same organism (same TaxID) share
the same key across different chain CSV files. This enables Boltz-2 to
perform evolutionary covariance analysis across interacting chains.
Pairing itself is not done by this script or anywhere else in this
pipeline: Boltz-2's own featurizer groups rows sharing the same key
across whichever two chains' CSVs are loaded for a given prediction,
purely by value equality.

TaxID lookup follows Boltz-2's own training-time rule exactly: a row only
participates in pairing if its header starts with "UniRef100" (after the
leading ">"); anything else (UniRef90, stripped/ID-only headers, templates,
etc.) gets key=-1 unconditionally. Unlike Boltz-2's own code (which reads an
embedded TaxID straight out of the header text), this script does an
EXTERNAL lookup against a running redis instance loaded with Boltz's
own taxonomy.rdb, keyed by the raw UniRef100 accession.

Sequences without a parseable accession, or whose accession isn't found
in the taxonomy lookup, receive key -1 (unpaired). The query sequence
(first entry) always receives key 0.

Input:  a3m file (from mmseqs2, colabfold_search, or jackhmmer)
Output: CSV file suitable for Boltz-2 YAML msa field
"""

import argparse
import csv
import sys

import redis


def parse_a3m(a3m_path):
    """
    Parse an a3m file, yielding (header, sequence) tuples.

    Lowercase letters (a3m convention for insertions) are preserved as-is,
    not stripped: Boltz-2's own CSV reader (boltz/data/parse/csv.py)
    reconstructs the deletion matrix by counting consecutive lowercase
    characters between aligned (uppercase/'-') columns -- identical logic to
    its native .a3m parser (boltz/data/parse/a3m.py). Stripping lowercase
    here would silently discard that deletion-matrix signal before it ever
    reaches Boltz-2, even though the CSV format itself fully supports it.
    """
    header = None
    seq_lines = []
    with open(a3m_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line[1:]
                seq_lines = []
            elif line.startswith("#"):
                continue  # skip comment lines
            elif line:
                seq_lines.append(line)
    if header is not None:
        yield header, "".join(seq_lines)


def extract_accession(header):
    """
    Extract the UniRef100 accession from a header, using Boltz-2's own
    eligibility rule: only headers starting with "UniRef100" are ever
    looked up; everything else gets key -1 unconditionally.

    The accession is the second underscore-delimited token; anything after
    it (e.g. mmseqs2's --msa-format-mode 6 stats, tab-separated) is
    discarded.
    """
    if not header.startswith("UniRef100"):
        return None
    parts = header.split("_")
    if len(parts) < 2:
        return None
    return parts[1].split("\t")[0].split()[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("a3m", help="Input a3m file")
    parser.add_argument("--output", "-o", required=True, help="Output CSV file")
    parser.add_argument("--max-seqs", type=int, default=8192,
                        help="Maximum number of MSA sequences to retain (default: 8192)")
    parser.add_argument("--redis-host", default="localhost",
                        help="Host of the redis instance loaded with taxonomy.rdb")
    parser.add_argument("--redis-port", type=int, default=7777,
                        help="Port of the redis instance loaded with taxonomy.rdb")
    args = parser.parse_args()

    r = redis.Redis(host=args.redis_host, port=args.redis_port)

    rows_written = 0
    n_eligible = 0
    n_found = 0

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
                acc = extract_accession(header)
                if acc is None:
                    key = -1
                else:
                    n_eligible += 1
                    val = r.get(acc)
                    if val is None:
                        key = -1
                    else:
                        n_found += 1
                        key = int(val)

            writer.writerow([seq, key])
            rows_written += 1

    print(f"Converted {rows_written} sequences from {args.a3m} -> {args.output} "
          f"({n_eligible} UniRef100-eligible, {n_found} found in taxonomy.rdb)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
