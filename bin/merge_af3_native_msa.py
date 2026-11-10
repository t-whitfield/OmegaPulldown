#!/usr/bin/env python3
"""
Merge per-database a3m MSA files (from the AF3-native unpaired-MSA search:
UniRef90, MGnify, Small BFD) into a single unpaired MSA per protein. UniProt
is deliberately NOT one of the inputs here: it is reserved for the
separate paired-MSA path (see MSA_SEARCH_UNIPROT_PAIRED).

Each input file was searched against a different database but shares the
same query sequence as its first record. The query record is kept once
(from the first input file); every subsequent record from every input file
is appended, with exact-duplicate sequences removed (ignoring lowercase
insertion characters). This matches AlphaFold 3's own native merge exactly.
"""

import argparse


def parse_a3m_records(path):
    """Yield (header_line, sequence_line) pairs, preserving original text."""
    header = None
    seq_lines = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "\n".join(seq_lines)
                header = line
                seq_lines = []
            elif line.startswith("#"):
                continue
            else:
                seq_lines.append(line)
    if header is not None:
        yield header, "\n".join(seq_lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("a3m_files", nargs="+", help="Per-database a3m files for one protein")
    parser.add_argument("--output", "-o", required=True, help="Merged output a3m file")
    args = parser.parse_args()

    query_written = False
    seen_sequences = set()
    with open(args.output, "w") as out:
        for i, path in enumerate(args.a3m_files):
            for j, (header, seq) in enumerate(parse_a3m_records(path)):
                is_query_record = (j == 0)
                if is_query_record:
                    if query_written:
                        continue
                    query_written = True
                else:
                    # Exact-duplicate dedup, ignoring lowercase insertions
                    # (a3m convention) -- matches AF3-native's own
                    # Msa(deduplicate=True) semantics.
                    dedup_key = "".join(ch for ch in seq if not ch.islower())
                    if dedup_key in seen_sequences:
                        continue
                    seen_sequences.add(dedup_key)
                out.write(f"{header}\n{seq}\n")


if __name__ == "__main__":
    main()
