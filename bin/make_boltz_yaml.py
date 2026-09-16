#!/usr/bin/env python3
"""
Generate a Boltz-2 YAML input file for a protein-protein pair.

Creates the YAML specification that Boltz-2 requires, optionally referencing
pre-computed MSA CSV files and post-translational modifications for each chain.

YAML format (with MSA and modifications):
  version: 1
  sequences:
    - protein:
        id: A
        sequence: <bait_sequence>
        msa: <path_to_bait_msa.csv>
        modifications:
          - position: 18
            ccd: SEP
    - protein:
        id: B
        sequence: <candidate_sequence>
        msa: <path_to_candidate_msa.csv>

Post-translational modifications are specified via optional CSV files
(one per role) with columns: protein_id, position, ccd

Common CCD codes:
  SEP  = phosphoserine        TPO = phosphothreonine
  PTR  = phosphotyrosine      MLY = N-dimethyl-lysine
  M3L  = N-trimethyl-lysine   ALY = N-acetyl-lysine
  HYP  = hydroxyproline       OCS = S-sulfo-cysteine
"""

import argparse
import csv
import json
import os
import sys
import yaml


def load_ptms_from_csv(csv_path, protein_id):
    """Load PTMs for a specific protein from a CSV file.

    CSV format: protein_id,position,ccd
    Returns a list of {"position": int, "ccd": str} dicts, or empty list.
    """
    if csv_path is None or csv_path.lower() in ("none", "null", "no_file", ""):
        return []

    if not os.path.isfile(csv_path):
        return []

    ptms = []
    with open(csv_path) as fh:
        reader = csv.DictReader(fh)

        # Validate header
        required = {"protein_id", "position", "ccd"}
        if not required.issubset(set(reader.fieldnames or [])):
            print(f"WARNING: PTM CSV {csv_path} missing required columns "
                  f"{required - set(reader.fieldnames or [])}; skipping",
                  file=sys.stderr)
            return []

        for row in reader:
            if row["protein_id"].strip() != protein_id:
                continue
            try:
                pos = int(row["position"])
            except ValueError:
                print(f"WARNING: non-integer position '{row['position']}' "
                      f"for {protein_id} in {csv_path}; skipping row",
                      file=sys.stderr)
                continue
            ccd = row["ccd"].strip()
            if not ccd:
                print(f"WARNING: empty CCD code at position {pos} "
                      f"for {protein_id} in {csv_path}; skipping row",
                      file=sys.stderr)
                continue
            ptms.append({"position": pos, "ccd": ccd})

    if ptms:
        print(f"  {protein_id}: {len(ptms)} PTM(s) loaded from {csv_path}",
              file=sys.stderr)
    return ptms


def validate_ptms(ptms, sequence, label):
    """Validate that PTM positions fall within the sequence."""
    seq_len = len(sequence)
    for ptm in ptms:
        pos = ptm["position"]
        if pos < 1 or pos > seq_len:
            print(f"WARNING: {label} PTM position {pos} (ccd={ptm['ccd']}) "
                  f"is outside sequence length {seq_len}",
                  file=sys.stderr)


def build_protein_entry(chain_id, sequence, msa_path=None, ptms=None):
    """Build a single protein entry for the Boltz-2 YAML."""
    entry = {
        "id": chain_id,
        "sequence": sequence,
    }
    if msa_path is not None:
        entry["msa"] = os.path.basename(msa_path)
    if ptms:
        entry["modifications"] = ptms
    return {"protein": entry}


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bait-id", required=True, help="Bait protein ID")
    parser.add_argument("--bait-seq", required=True, help="Bait protein sequence")
    parser.add_argument("--bait-msa", default=None, help="Path to bait MSA CSV")
    parser.add_argument("--candidate-id", required=True, help="Candidate protein ID")
    parser.add_argument("--candidate-seq", required=True, help="Candidate protein sequence")
    parser.add_argument("--candidate-msa", default=None, help="Path to candidate MSA CSV")
    parser.add_argument("--bait-ptms", default=None,
                        help="Path to bait PTM CSV (columns: protein_id,position,ccd)")
    parser.add_argument("--candidate-ptms", default=None,
                        help="Path to candidate PTM CSV (columns: protein_id,position,ccd)")
    parser.add_argument("--output", "-o", required=True, help="Output YAML file")
    args = parser.parse_args()

    # Validate sequences contain only standard amino acids
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    for label, seq in [("bait", args.bait_seq), ("candidate", args.candidate_seq)]:
        invalid = set(seq.upper()) - valid_aa
        if invalid:
            print(f"WARNING: {label} sequence contains non-standard residues: {invalid}",
                  file=sys.stderr)

    # Load PTMs from optional CSV files
    bait_ptms = load_ptms_from_csv(args.bait_ptms, args.bait_id)
    candidate_ptms = load_ptms_from_csv(args.candidate_ptms, args.candidate_id)

    # Validate PTM positions against sequences
    validate_ptms(bait_ptms, args.bait_seq, "bait")
    validate_ptms(candidate_ptms, args.candidate_seq, "candidate")

    boltz_input = {
        "version": 1,
        "sequences": [
            build_protein_entry("A", args.bait_seq, args.bait_msa, bait_ptms),
            build_protein_entry("B", args.candidate_seq, args.candidate_msa, candidate_ptms),
        ]
    }

    with open(args.output, "w") as fh:
        yaml.dump(boltz_input, fh, default_flow_style=False, sort_keys=False)

    print(f"Wrote Boltz-2 YAML: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
