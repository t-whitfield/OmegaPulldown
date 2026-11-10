#!/usr/bin/env python3
"""
pDockQ from a predicted complex structure.

pDockQ follows Bryant et al. (2022, Nat Commun 13:6028):
    x      = avg_interface_pLDDT * log10(n_interface_contacts)
    pDockQ = L / (1 + exp(-k * (x - x0))) + b
with L=0.724, x0=152.611, k=0.052, b=0.018.

Interface contacts are Cbeta-Cbeta pairs below 10 A between the first two
chains (Calpha standing in for Cbeta at glycine). avg_interface_pLDDT is the
mean over the union of interface residues on both chains: residue-weighted,
not contact-weighted.

Importable:
    from compute_pdockq import compute_pdockq
    pdockq, n_contacts, avg_if_plddt = compute_pdockq("model.cif")
    # pdockq is None when the chains make no contact, or there are <2 chains

Runnable, for spot-checks outside the pipeline:
    compute_pdockq.py model_0.cif model_1.cif ...
    -> TSV: cif, pdockq, n_contacts, avg_interface_plddt  (pdockq empty if None)
"""

import math
import sys

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
        # None, not 0.0 -- "no interface" is not the same measurement as "an
        # interface that scores zero", and the callers rely on the
        # distinction: PARSE_METRICS skips None rather than averaging it into
        # pdockq_mean. This matches the behaviour the embedded copies have
        # always had; the standalone script this was promoted from returned
        # 0.0 here, and that difference was reconciled in favour of None.
        return None, 0, 0.0

    # Interface residue indices
    a_interface = np.any(contact_mask, axis=1)
    b_interface = np.any(contact_mask, axis=0)

    interface_plddts = np.concatenate([plddts_a[a_interface], plddts_b[b_interface]])
    avg_interface_plddt = float(np.mean(interface_plddts))

    # pDockQ formula
    x = avg_interface_plddt * math.log10(n_contacts)
    pdockq = PDOCKQ_L / (1.0 + math.exp(-PDOCKQ_K * (x - PDOCKQ_X0))) + PDOCKQ_B

    return pdockq, n_contacts, avg_interface_plddt


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__.strip().splitlines()[0])
    print("cif\tpdockq\tn_contacts\tavg_interface_plddt")
    for path in sys.argv[1:]:
        pdockq, n_contacts, avg_if_plddt = compute_pdockq(path)
        print("%s\t%s\t%d\t%.4f" % (
            path, "" if pdockq is None else "%.6f" % pdockq,
            n_contacts, avg_if_plddt))


if __name__ == "__main__":
    main()
