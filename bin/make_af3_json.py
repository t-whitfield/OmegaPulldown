#!/usr/bin/env python3
"""
Generate an AlphaFold 3 JSON input file for a protein-protein pair.

Creates a version-2 JSON that references pre-computed A3M MSA files via
unpairedMsaPath, allowing AF3 to run with --norun_data_pipeline.

Each protein chain specifies:
  - unpairedMsaPath: relative path to the A3M file
  - pairedMsa: "" (empty; no paired MSA)
  - templates: [] (empty; no structural templates)

Post-translational modifications use ptmType/ptmPosition fields
in the AF3 JSON format.

Usage:
  make_af3_json.py \
      --bait-id Wnt3a --bait-seq SYPIWW... --bait-msa Wnt3a.a3m \
      --candidate-id xCres --candidate-seq MYLDFF... --candidate-msa xCres.a3m \
      --bait-ptms bait_ptms.csv \
      --seeds 42 --output Wnt3a_with_xCres.json
"""

import argparse
import csv
import json
import os
import sys


def load_ptms_from_csv(csv_path, protein_id):
    """Load PTMs for a specific protein from a CSV file.

    CSV format: protein_id,position,ccd
    Returns a list of {"ptmType": str, "ptmPosition": int} dicts for AF3.
    """
    if csv_path is None or csv_path.lower() in ("none", "null", "no_file", ""):
        return []

    if not os.path.isfile(csv_path):
        return []

    ptms = []
    with open(csv_path) as fh:
        reader = csv.DictReader(fh)

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
            ptms.append({"ptmType": ccd, "ptmPosition": pos})

    if ptms:
        print(f"  {protein_id}: {len(ptms)} PTM(s) loaded from {csv_path}",
              file=sys.stderr)
    return ptms


def validate_ptms(ptms, sequence, label):
    """Validate that PTM positions fall within the sequence."""
    seq_len = len(sequence)
    for ptm in ptms:
        pos = ptm["ptmPosition"]
        if pos < 1 or pos > seq_len:
            print(f"WARNING: {label} PTM position {pos} (ccd={ptm['ptmType']}) "
                  f"is outside sequence length {seq_len}",
                  file=sys.stderr)


def clean_a3m_for_af3(a3m_path, output_path):
    """Clean an A3M file for AF3 compatibility.

    Removes characters that AF3's MSA featuriser rejects:
      - X (ambiguous residue) in hit sequences → replaced with gap (-)
      - Null bytes and other control characters → stripped
      - Trailing blank/corrupt lines → removed

    The query (first) sequence is NOT modified here; use patch_a3m_query()
    for that.

    Returns the number of substitutions made.
    """
    with open(a3m_path) as f:
        lines = f.readlines()

    cleaned = []
    query_done = False
    in_query = False
    n_subs = 0
    for line in lines:
        # Strip null bytes and control chars (except newline)
        line = "".join(ch for ch in line if ch == "\n" or (ch >= " " and ch != "\x7f"))
        if not line.strip():
            continue  # skip empty lines

        if line.startswith(">"):
            if not in_query and not query_done:
                in_query = True
            elif in_query:
                in_query = False
                query_done = True
            cleaned.append(line)
        elif in_query:
            # Query sequence: pass through unmodified
            cleaned.append(line)
        else:
            # Hit sequence: replace X with gap
            new_line = line.replace("X", "-")
            n_subs += line.count("X")
            cleaned.append(new_line)

    with open(output_path, "w") as out:
        out.writelines(cleaned)

    if n_subs > 0:
        print(f"  Cleaned A3M for AF3: replaced {n_subs} X chars in hits → {output_path}",
              file=sys.stderr)
    return n_subs


def patch_a3m_query(a3m_path, ptm_positions, output_path):
    """Rewrite the A3M file, replacing modified positions with X in the query.

    AF3 internally substitutes modified residues with X in its query sequence
    and requires the MSA's first (query) sequence to match.  The A3M from
    mmseqs2 (or jackhmmer) has the original unmodified sequence, so we patch
    it here.

    Parameters
    ----------
    a3m_path : str
        Path to the original A3M file.
    ptm_positions : set of int
        1-based positions to replace with X.
    output_path : str
        Path to write the patched A3M file.
    """
    with open(a3m_path) as fh:
        lines = fh.readlines()

    if not ptm_positions:
        # No modifications — just copy the file
        with open(output_path, "w") as out:
            out.writelines(lines)
        return

    # The query sequence is the first sequence after the first header line.
    # It may span multiple lines before the next header (line starting with >).
    patched = []
    query_done = False
    in_query = False
    for line in lines:
        if line.startswith(">") and not in_query and not query_done:
            patched.append(line)
            in_query = True
            seq_pos = 0  # track 1-based position in the alignment (uppercase only)
            continue
        if in_query:
            if line.startswith(">"):
                # End of query sequence, start of next entry
                in_query = False
                query_done = True
                patched.append(line)
                continue
            # Patch this line of the query sequence
            new_chars = []
            for ch in line.rstrip("\n"):
                if ch.isupper() or ch == "-":
                    seq_pos += 1
                    if ch.isupper() and seq_pos in ptm_positions:
                        new_chars.append("X")
                    else:
                        new_chars.append(ch)
                else:
                    # Lowercase = insertion; don't increment position
                    new_chars.append(ch)
            patched.append("".join(new_chars) + "\n")
        else:
            patched.append(line)

    with open(output_path, "w") as out:
        out.writelines(patched)

    print(f"  Patched A3M query: {len(ptm_positions)} position(s) -> X in {output_path}",
          file=sys.stderr)


def build_protein_entry(chain_id, sequence, msa_path=None, ptms=None):
    """Build a single protein entry for the AF3 JSON."""
    entry = {
        "id": chain_id,
        "sequence": sequence,
    }
    if ptms:
        entry["modifications"] = ptms
    if msa_path is not None:
        entry["unpairedMsaPath"] = msa_path
        entry["pairedMsa"] = ""
        entry["templates"] = []
    return {"protein": entry}


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bait-id", required=True, help="Bait protein ID")
    parser.add_argument("--bait-seq", required=True, help="Bait protein sequence")
    parser.add_argument("--bait-msa", default=None, help="Path to bait A3M MSA file")
    parser.add_argument("--candidate-id", required=True, help="Candidate protein ID")
    parser.add_argument("--candidate-seq", required=True, help="Candidate protein sequence")
    parser.add_argument("--candidate-msa", default=None, help="Path to candidate A3M MSA file")
    parser.add_argument("--bait-ptms", default=None,
                        help="Path to bait PTM CSV (columns: protein_id,position,ccd)")
    parser.add_argument("--candidate-ptms", default=None,
                        help="Path to candidate PTM CSV (columns: protein_id,position,ccd)")
    parser.add_argument("--seeds", default="42",
                        help="Comma-separated model seeds (default: 42)")
    parser.add_argument("--output", "-o", required=True, help="Output JSON file")
    args = parser.parse_args()

    # Validate sequences
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    for label, seq in [("bait", args.bait_seq), ("candidate", args.candidate_seq)]:
        invalid = set(seq.upper()) - valid_aa
        if invalid:
            print(f"WARNING: {label} sequence contains non-standard residues: {invalid}",
                  file=sys.stderr)

    # Parse model seeds
    seeds = [int(s.strip()) for s in args.seeds.split(",")]

    # Load PTMs
    bait_ptms = load_ptms_from_csv(args.bait_ptms, args.bait_id)
    candidate_ptms = load_ptms_from_csv(args.candidate_ptms, args.candidate_id)

    validate_ptms(bait_ptms, args.bait_seq, "bait")
    validate_ptms(candidate_ptms, args.candidate_seq, "candidate")

    # --- Patch A3M query sequences at PTM positions (AF3 expects X) ---
    # AF3 internally replaces modified residues with X in its query sequence
    # and requires the MSA's first (query) sequence to match.  We patch the
    # A3M files here so they agree with AF3's internal representation.

    # Determine MSA references for each chain.
    # For chains WITH PTMs: embed the MSA content inline via "unpairedMsa"
    # (a string field) with the query sequence patched to use X at modified
    # positions.  AF3's JSON parser requires the query to match its internal
    # X-substituted sequence, but the file-based "unpairedMsaPath" code path
    # also feeds X into the MSA featuriser, which rejects it.  The inline
    # "unpairedMsa" path avoids this by letting AF3 re-derive features from
    # the raw alignment string while still passing the query-sequence check.
    #
    # For chains WITHOUT PTMs: use "unpairedMsaPath" as before.

    # --- Clean and patch A3M files for AF3 ---
    # Two issues with mmseqs2/jackhmmer A3Ms that AF3 rejects:
    #   1. Hit sequences may contain X (ambiguous residue) → replace with gap
    #   2. Files may have null bytes or trailing corruption → strip
    # Additionally, for PTM chains, AF3 requires the query sequence to have
    # X at modified positions (AF3 does this internally and checks the MSA).
    #
    # Strategy: clean ALL A3Ms, then patch query for PTM chains.

    bait_msa_ref = os.path.basename(args.bait_msa) if args.bait_msa else None
    if args.bait_msa:
        cleaned_name = os.path.basename(args.bait_msa).replace(".a3m", "_af3.a3m")
        clean_a3m_for_af3(args.bait_msa, cleaned_name)
        if bait_ptms:
            ptm_positions = {p["ptmPosition"] for p in bait_ptms}
            patch_a3m_query(cleaned_name, ptm_positions, cleaned_name)
        bait_msa_ref = cleaned_name

    cand_msa_ref = os.path.basename(args.candidate_msa) if args.candidate_msa else None
    if args.candidate_msa:
        cleaned_name = os.path.basename(args.candidate_msa).replace(".a3m", "_af3.a3m")
        clean_a3m_for_af3(args.candidate_msa, cleaned_name)
        if candidate_ptms:
            ptm_positions = {p["ptmPosition"] for p in candidate_ptms}
            patch_a3m_query(cleaned_name, ptm_positions, cleaned_name)
        cand_msa_ref = cleaned_name

    # Build AF3 JSON (version 2 enables external MSA fields)
    af3_input = {
        "dialect": "alphafold3",
        "version": 2,
        "name": args.output.replace(".json", ""),
        "modelSeeds": seeds,
        "sequences": [
            build_protein_entry(
                "A", args.bait_seq,
                msa_path=bait_msa_ref,
                ptms=bait_ptms if bait_ptms else None,
            ),
            build_protein_entry(
                "B", args.candidate_seq,
                msa_path=cand_msa_ref,
                ptms=candidate_ptms if candidate_ptms else None,
            ),
        ],
    }

    with open(args.output, "w") as fh:
        json.dump(af3_input, fh, indent=2)

    print(f"Wrote AF3 JSON: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
