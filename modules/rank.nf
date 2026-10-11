/*
 * Parse Boltz-2 outputs and rank predictions by interface quality.
 *
 * Two-step approach:
 *   1. PARSE_METRICS: per-pair process extracting one row of metrics
 *   2. RANK_PREDICTIONS: collects all rows, sorts, writes final TSV
 *
 * This is more Nextflow-idiomatic than collecting directories.
 */

process PARSE_METRICS {
    tag "${pair_id}"
    label 'process_low'

    input:
    tuple val(pair_id), path(pred_files)

    output:
    path "${pair_id}.metrics.tsv", emit: metrics

    script:
    """
    #!/usr/bin/env python3
    import json, glob, math, re, sys, os
    import numpy as np

    pair_id = "${pair_id}"

    # =====================================================================
    # Helper: parse one CIF and return pDockQ, n_contacts, avg_if_plddt
    # =====================================================================
    def compute_pdockq(cif_path):
        atoms = []
        in_atom_site = False
        columns = []
        with open(cif_path) as f:
            for line in f:
                line = line.strip()
                if line.startswith("_atom_site."):
                    in_atom_site = True
                    columns.append(line.split(".")[1].strip())
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
                    rec = dict(zip(columns, tokens))
                    if rec.get("group_PDB") != "ATOM":
                        continue
                    try:
                        atoms.append({
                            "atom": rec.get("label_atom_id",""),
                            "chain": rec.get("label_asym_id",""),
                            "resi": int(rec.get("label_seq_id","0")),
                            "x": float(rec.get("Cartn_x","0")),
                            "y": float(rec.get("Cartn_y","0")),
                            "z": float(rec.get("Cartn_z","0")),
                            "b": float(rec.get("B_iso_or_equiv","0")),
                        })
                    except (ValueError, KeyError):
                        pass

        ca, cb = {}, {}
        for a in atoms:
            k = (a["chain"], a["resi"])
            if a["atom"] == "CA": ca[k] = (a["x"],a["y"],a["z"],a["b"])
            elif a["atom"] == "CB": cb[k] = (a["x"],a["y"],a["z"],a["b"])
        rep = {k: cb.get(k, ca[k]) for k in ca}

        chains = sorted(set(k[0] for k in rep))
        if len(chains) < 2:
            return None, 0, 0.0
        ra = {k:v for k,v in rep.items() if k[0]==chains[0]}
        rb = {k:v for k,v in rep.items() if k[0]==chains[1]}
        if not ra or not rb:
            return None, 0, 0.0
        ca_ = np.array([[v[0],v[1],v[2]] for v in ra.values()])
        cb_ = np.array([[v[0],v[1],v[2]] for v in rb.values()])
        pa = np.array([v[3] for v in ra.values()])
        pb = np.array([v[3] for v in rb.values()])
        d = np.sqrt(np.sum((ca_[:,None,:]-cb_[None,:,:])**2, axis=2))
        mask = d < 10.0
        n_contacts = int(np.sum(mask))
        if n_contacts == 0:
            return None, 0, 0.0
        ai = np.any(mask, axis=1)
        bi = np.any(mask, axis=0)
        avg_if_plddt = float(np.mean(np.concatenate([pa[ai],pb[bi]])))
        x = avg_if_plddt * math.log10(n_contacts)
        pdockq = 0.724/(1+math.exp(-0.052*(x-152.611)))+0.018
        return pdockq, n_contacts, avg_if_plddt

    # =====================================================================
    # Collect per-sample metrics
    # =====================================================================
    conf_keys = ['confidence_score','ptm','iptm','protein_iptm',
                 'complex_plddt','complex_iplddt']
    struct_keys = ['pdockq','n_interface_contacts','avg_interface_plddt']
    all_keys = conf_keys + struct_keys

    # Discover samples — Boltz-2 names them *_model_0.json .. *_model_N.json
    conf_files = sorted(glob.glob("confidence_*.json"))
    cif_files  = sorted(glob.glob("*.cif"))
    n_samples  = max(len(conf_files), len(cif_files), 1)

    # Accumulate per-sample values for each metric
    sample_vals = {k: [] for k in all_keys}

    for i in range(n_samples):
        # --- Confidence metrics ---
        if i < len(conf_files):
            with open(conf_files[i]) as f:
                data = json.load(f)
            for key in conf_keys:
                v = data.get(key)
                if v is not None:
                    sample_vals[key].append(float(v))
        # --- Structural metrics (pDockQ) ---
        if i < len(cif_files):
            pdockq, nc, aip = compute_pdockq(cif_files[i])
            if pdockq is not None:
                sample_vals['pdockq'].append(pdockq)
            sample_vals['n_interface_contacts'].append(float(nc))
            sample_vals['avg_interface_plddt'].append(aip)

    if not conf_files and not cif_files:
        print(f"WARNING: no confidence JSON or CIF for {pair_id}", file=sys.stderr)

    # =====================================================================
    # Compute mean and std, write TSV
    # =====================================================================
    out_cols = ['pair_id', 'n_samples']
    for k in all_keys:
        out_cols.append(f"{k}_mean")
        out_cols.append(f"{k}_std")

    out_vals = [pair_id, str(n_samples)]
    for k in all_keys:
        arr = sample_vals[k]
        if len(arr) == 0:
            out_vals.extend(["NA", "NA"])
        else:
            a = np.array(arr)
            out_vals.append(f"{np.mean(a):.4f}")
            out_vals.append(f"{np.std(a, ddof=0):.4f}" if len(a) > 1 else "NA")

    with open("${pair_id}.metrics.tsv", "w") as out:
        out.write("\\t".join(out_cols) + "\\n")
        out.write("\\t".join(out_vals) + "\\n")

    print(f"{pair_id}: {n_samples} samples, "
          f"protein_iptm = {np.mean(sample_vals['protein_iptm']):.4f} "
          f"+/- {np.std(sample_vals['protein_iptm'], ddof=0):.4f}"
          if sample_vals['protein_iptm'] else f"{pair_id}: no confidence data",
          file=sys.stderr)
    """
}


process RANK_PREDICTIONS {
    label 'process_low'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path "*"   // all per-pair metric TSV files

    output:
    path "rankings_${params.predictor}.tsv", emit: rankings

    script:
    """
    #!/usr/bin/env python3
    import glob, sys

    predictor = "${params.predictor}"
    header = None
    rows = []

    for tsv in sorted(glob.glob("*.metrics.tsv")):
        with open(tsv) as f:
            h = f.readline().strip()
            if header is None:
                header = h
            row = f.readline().strip()
            if row:
                rows.append(row.split("\\t"))

    # The rank_by param names the base metric (e.g. protein_iptm).
    # Columns are now named {metric}_mean, {metric}_std, so append _mean.
    rank_by = "${params.rank_by}"
    cols = header.split("\\t")
    rank_col = rank_by + "_mean"
    try:
        idx = cols.index(rank_col)
    except ValueError:
        # Fallback: try the raw name, then default to column 2 (first _mean col)
        try:
            idx = cols.index(rank_by)
        except ValueError:
            idx = 2

    def sortkey(r):
        try:
            return float(r[idx])
        except (ValueError, IndexError):
            return -1.0

    rows.sort(key=sortkey, reverse=True)

    outname = f"rankings_{predictor}.tsv"
    with open(outname, "w") as out:
        out.write("rank\\t" + header + "\\n")
        for i, r in enumerate(rows, 1):
            out.write(f"{i}\\t" + "\\t".join(r) + "\\n")

    print(f"Ranked {len(rows)} predictions by {rank_col}", file=sys.stderr)
    if rows:
        print(f"Top: {rows[0][0]} ({rank_col}={sortkey(rows[0]):.4f})", file=sys.stderr)
    """
}
