#!/usr/bin/env python3
"""
Convert mmseqs2 result2msa --msa-format-mode 2 ("aligned FASTA DB") output
into A3M format, preserving full original database headers.

Mode 2 is used (instead of mode 6, used elsewhere in OmegaPulldown) only for
the UniProt paired-MSA search: mode 6 discards the original header text in
favor of a synthetic ID+stats header, but AlphaFold 3's own species-pairing
code needs the raw UniProt header (`sp|ACCESSION|GENE_SPECIES` /
`tr|ACCESSION|GENE_SPECIES`) to extract a per-hit species identifier.
"""

import argparse


def parse_fasta_records(path):
    """Yield (header_without_gt, sequence) pairs from a FASTA-like file."""
    header = None
    seq_lines = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        yield header, "".join(seq_lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("aligned_fasta", help="mmseqs2 --msa-format-mode 2 output")
    parser.add_argument("--output", "-o", required=True, help="Output A3M file")
    args = parser.parse_args()

    with open(args.output, "w") as out:
        for i, (header, seq) in enumerate(parse_fasta_records(args.aligned_fasta)):
            if i == 0:
                seq = seq.replace("-", "")
            out.write(f">{header}\n{seq}\n")


if __name__ == "__main__":
    main()
