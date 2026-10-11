#!/usr/bin/env python3
"""
Parse Boltz-2 prediction outputs and compute ranking metrics.

For each predicted complex, extracts:
  - Boltz-2 native confidence metrics from confidence JSON
    (ipTM, protein_ipTM, pTM, complex_pLDDT, complex_ipLDDT, confidence_score,
     pair_chains_iptm)
  - pDockQ score computed from the predicted CIF structure

pDockQ calculation follows Bryant et al. (2022, Nat Commun 13:6028):
  x = avg_interface_pLDDT * log10(n_interface_contacts)
  pDockQ = L / (1 + exp(-k*(x - x0))) + b
  where L=0.724, x0=152.611, k=0.052, b=0.018

Interface contacts: Cbeta-Cbeta distance < 10 Angstrom between chains
(Calpha for glycine).

Output: ranked TSV with all metrics, sorted by chosen ranking metric.
"""

import argparse
import json
import math
import os
import sys
from pathlib import Path

import numpy as np

# pDockQ sigmoid parameters (Bryant et al. 2022)
PDOCKQ_L = 0.724
PDOCKQ_X0 = 152.611
PDOCKQ_K = 0.052
PDOCKQ_B = 0.018

# Interface contact distance threshold (Angstrom)
INTERFACE_CUTOFF = 10.0


def parse_cif_atoms(cif_path):
    """
    Minimal CIF parser: extract atom coordinates, chain IDs, residue info,
    and B-factor (pLDDT) from _atom_site records.

    Returns list of dicts with keys:
      atom_name, chain_id, res_name, res_seq, x, y, z, bfactor
    """
    atoms = []
    in_atom_site = False
    columns = []

    with open(cif_path) as fh:
        for line in fh:
            line = line.strip()

            if line.startswith("_atom_site."):
                in_atom_site = True
                col_name = line.split(".")[1].strip()
                columns.append(col_name)
                continue

            if in_atom_site and (line.startswith("_") or line.startswith("#") or line.startswith("loop_")):
                in_atom_site = False
                if line.startswith("loop_"):
                    columns = []
                continue

            if in_atom_site and columns and not line.startswith("_"):
                tokens = line.split()
                if len(tokens) < len(columns):
                    continue

                record = dict(zip(columns, tokens))
                # Only keep ATOM records (not HETATM unless part of protein)
                group = record.get("group_PDB", "ATOM")
                if group not in ("ATOM",):
                    continue

                try:
                    atom = {
                        "atom_name": record.get("label_atom_id", record.get("auth_atom_id", "")),
                        "chain_id": record.get("label_asym_id", record.get("auth_asym_id", "")),
                        "res_name": record.get("label_comp_id", record.get("auth_comp_id", "")),
                        "res_seq": int(record.get("label_seq_id", record.get("auth_seq_id", 0))),
                        "x": float(record.get("Cartn_x", 0)),
                        "y": float(record.get("Cartn_y", 0)),
                        "z": float(record.get("Cartn_z", 0)),
                        "bfactor": float(record.get("B_iso_or_equiv", 0)),
                    }
                    atoms.append(atom)
                except (ValueError, KeyError):
                    continue

    return atoms


def get_cb_atoms(atoms):
    """
    Extract Cbeta atoms (or Calpha for glycine) per residue per chain.
    Returns dict: {(chain_id, res_seq): (x, y, z, plddt)}
    """
    # First pass: collect all Calpha and Cbeta
    ca_atoms = {}
    cb_atoms = {}
    for atom in atoms:
        key = (atom["chain_id"], atom["res_seq"])
        if atom["atom_name"] == "CA":
            ca_atoms[key] = (atom["x"], atom["y"], atom["z"], atom["bfactor"])
        elif atom["atom_name"] == "CB":
            cb_atoms[key] = (atom["x"], atom["y"], atom["z"], atom["bfactor"])

    # Use CB where available, CA (glycine) otherwise
    result = {}
    for key in ca_atoms:
        if key in cb_atoms:
            result[key] = cb_atoms[key]
        else:
            result[key] = ca_atoms[key]
    return result


def compute_pdockq(cif_path):
    """
    Compute pDockQ from a predicted CIF file.

    Returns (pdockq, n_contacts, avg_interface_plddt) or (None, 0, 0)
    if computation fails.
    """
    atoms = parse_cif_atoms(cif_path)
    if not atoms:
        return None, 0, 0.0

    cb = get_cb_atoms(atoms)

    # Separate chains
    chains = sorted(set(k[0] for k in cb.keys()))
    if len(chains) < 2:
        return None, 0, 0.0

    # For binary complex: chain A vs chain B
    chain_a = chains[0]
    chain_b = chains[1]

    residues_a = {k: v for k, v in cb.items() if k[0] == chain_a}
    residues_b = {k: v for k, v in cb.items() if k[0] == chain_b}

    # Compute interface contacts
    interface_plddts = []
    n_contacts = 0

    coords_a = np.array([[v[0], v[1], v[2]] for v in residues_a.values()])
    coords_b = np.array([[v[0], v[1], v[2]] for v in residues_b.values()])
    plddts_a = np.array([v[3] for v in residues_a.values()])
    plddts_b = np.array([v[3] for v in residues_b.values()])

    if len(coords_a) == 0 or len(coords_b) == 0:
        return None, 0, 0.0

    # Distance matrix: |A| x |B|
    diff = coords_a[:, np.newaxis, :] - coords_b[np.newaxis, :, :]
    dist_matrix = np.sqrt(np.sum(diff**2, axis=2))

    # Find contacts below threshold
    contact_mask = dist_matrix < INTERFACE_CUTOFF
    n_contacts = int(np.sum(contact_mask))

    if n_contacts == 0:
        return 0.0, 0, 0.0

    # Interface residue indices
    a_interface = np.any(contact_mask, axis=1)
    b_interface = np.any(contact_mask, axis=0)

    interface_plddts = np.concatenate([plddts_a[a_interface], plddts_b[b_interface]])
    avg_interface_plddt = float(np.mean(interface_plddts))

    # pDockQ formula
    x = avg_interface_plddt * math.log10(n_contacts)
    pdockq = PDOCKQ_L / (1.0 + math.exp(-PDOCKQ_K * (x - PDOCKQ_X0))) + PDOCKQ_B

    return pdockq, n_contacts, avg_interface_plddt


def parse_confidence_json(json_path):
    """Parse Boltz-2 confidence JSON and return metrics dict."""
    with open(json_path) as fh:
        data = json.load(fh)

    # Boltz-2 confidence JSON fields
    metrics = {
        "confidence_score": data.get("confidence_score", None),
        "ptm": data.get("ptm", None),
        "iptm": data.get("iptm", None),
        "protein_iptm": data.get("protein_iptm", None),
        "complex_plddt": data.get("complex_plddt", None),
        "complex_iplddt": data.get("complex_iplddt", None),
        "pair_chains_iptm": data.get("pair_chains_iptm", None),
    }

    return metrics


def process_prediction(pred_dir, pair_id):
    """
    Process a single Boltz-2 prediction directory.
    Returns a dict of all metrics for the best model.
    """
    pred_dir = Path(pred_dir)

    # Find confidence JSON (model 0 = best ranked)
    conf_files = sorted(pred_dir.glob("confidence_*_model_0.json"))
    if not conf_files:
        # Try alternative naming
        conf_files = sorted(pred_dir.glob("confidence_*.json"))
    if not conf_files:
        print(f"WARNING: no confidence JSON in {pred_dir}", file=sys.stderr)
        return None

    metrics = parse_confidence_json(conf_files[0])
    metrics["pair_id"] = pair_id

    # Find CIF structure
    cif_files = sorted(pred_dir.glob("*_model_0.cif"))
    if not cif_files:
        cif_files = sorted(pred_dir.glob("*.cif"))

    if cif_files:
        pdockq, n_contacts, avg_if_plddt = compute_pdockq(cif_files[0])
        metrics["pdockq"] = pdockq
        metrics["n_interface_contacts"] = n_contacts
        metrics["avg_interface_plddt"] = avg_if_plddt
    else:
        print(f"WARNING: no CIF structure in {pred_dir}", file=sys.stderr)
        metrics["pdockq"] = None
        metrics["n_interface_contacts"] = 0
        metrics["avg_interface_plddt"] = 0.0

    return metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--predictions-dir", required=True,
                        help="Directory containing Boltz-2 prediction outputs")
    parser.add_argument("--output", "-o", default="rankings.tsv",
                        help="Output TSV file with ranked results")
    parser.add_argument("--rank-by", default="protein_iptm",
                        choices=["protein_iptm", "iptm", "pdockq",
                                 "confidence_score", "complex_plddt"],
                        help="Metric to rank by (default: protein_iptm)")
    args = parser.parse_args()

    pred_base = Path(args.predictions_dir)

    # Collect results from all prediction subdirectories
    all_results = []

    # Boltz-2 output structure: predictions/<input_name>/
    pred_root = pred_base / "predictions" if (pred_base / "predictions").exists() else pred_base

    for subdir in sorted(pred_root.iterdir()):
        if not subdir.is_dir():
            continue
        pair_id = subdir.name
        result = process_prediction(subdir, pair_id)
        if result is not None:
            all_results.append(result)

    if not all_results:
        sys.exit("ERROR: no valid predictions found")

    # Sort by ranking metric (descending; higher is better)
    def sort_key(r):
        val = r.get(args.rank_by)
        return val if val is not None else -1.0

    all_results.sort(key=sort_key, reverse=True)

    # Write ranked output
    columns = [
        "rank", "pair_id", "protein_iptm", "iptm", "ptm",
        "confidence_score", "complex_plddt", "complex_iplddt",
        "pdockq", "n_interface_contacts", "avg_interface_plddt"
    ]

    with open(args.output, "w") as out:
        out.write("\t".join(columns) + "\n")
        for i, result in enumerate(all_results, 1):
            values = [str(i), result["pair_id"]]
            for col in columns[2:]:
                val = result.get(col)
                if val is None:
                    values.append("NA")
                elif isinstance(val, float):
                    values.append(f"{val:.4f}")
                else:
                    values.append(str(val))
            out.write("\t".join(values) + "\n")

    print(f"Ranked {len(all_results)} predictions -> {args.output}", file=sys.stderr)
    print(f"Top hit: {all_results[0]['pair_id']} "
          f"({args.rank_by}={sort_key(all_results[0]):.4f})", file=sys.stderr)


if __name__ == "__main__":
    main()
