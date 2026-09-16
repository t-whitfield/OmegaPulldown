/*
 * Parse AlphaFold 3 outputs and rank predictions by interface quality.
 *
 * AF3 output structure per pair:
 *   <pair_id>/
 *     <pair_id>_summary_confidences.json   (top-ranked model summary)
 *     <pair_id>_ranking_scores.csv         (per-sample ranking scores)
 *     <pair_id>_model.cif                  (top-ranked structure)
 *     seed-<N>_sample-<M>/
 *       <pair_id>_seed-<N>_sample-<M>_model.cif
 *       <pair_id>_seed-<N>_sample-<M>_summary_confidences.json
 *
 * Reuses RANK_PREDICTIONS from rank.nf — the output TSV format is identical.
 */

process PARSE_METRICS_AF3 {
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
    # Collect per-sample metrics from AF3 output
    # =====================================================================
    # AF3 metrics mapping to Boltz-2 column names for unified ranking:
    #   AF3 ranking_score -> confidence_score
    #   AF3 ptm           -> ptm
    #   AF3 iptm          -> iptm  (also used as protein_iptm)
    #   AF3 complex_plddt  -> from chain_ptm average
    #   AF3 has_clash      -> (logged but not ranked)

    conf_keys = ['confidence_score','ptm','iptm','protein_iptm',
                 'complex_plddt','complex_iplddt']
    struct_keys = ['pdockq','n_interface_contacts','avg_interface_plddt']
    all_keys = conf_keys + struct_keys

    sample_vals = {k: [] for k in all_keys}

    # Find per-sample confidence JSONs
    sample_conf_files = sorted(glob.glob("*_sample-*_summary_confidences.json"))
    sample_cif_files  = sorted(glob.glob("*_sample-*_model.cif"))

    # If no per-sample files, fall back to the top-level summary
    if not sample_conf_files:
        sample_conf_files = sorted(glob.glob("*_summary_confidences.json"))
        # Exclude any that contain "sample-" (shouldn't happen, but be safe)
        sample_conf_files = [f for f in sample_conf_files if "sample-" not in f]
    if not sample_cif_files:
        sample_cif_files = sorted(glob.glob("*_model.cif"))
        sample_cif_files = [f for f in sample_cif_files if "sample-" not in f]

    n_samples = max(len(sample_conf_files), len(sample_cif_files), 1)

    for i in range(n_samples):
        # --- Confidence metrics ---
        if i < len(sample_conf_files):
            with open(sample_conf_files[i]) as f:
                data = json.load(f)

            # AF3 -> unified column mapping
            if 'ranking_score' in data:
                sample_vals['confidence_score'].append(float(data['ranking_score']))
            if 'ptm' in data:
                sample_vals['ptm'].append(float(data['ptm']))
            if 'iptm' in data:
                sample_vals['iptm'].append(float(data['iptm']))
                # AF3 iptm is the inter-chain metric; map to protein_iptm for ranking
                sample_vals['protein_iptm'].append(float(data['iptm']))

            # complex_plddt: AF3 doesn't output this directly; approximate
            # from chain_ptm if available, otherwise leave as NA
            if 'chain_ptm' in data:
                avg_ptm = float(np.mean(data['chain_ptm']))
                sample_vals['complex_plddt'].append(avg_ptm)

            # complex_iplddt: not directly available in AF3
            # Leave empty (will be NA)

        # --- Structural metrics (pDockQ) ---
        if i < len(sample_cif_files):
            pdockq, nc, aip = compute_pdockq(sample_cif_files[i])
            if pdockq is not None:
                sample_vals['pdockq'].append(pdockq)
            sample_vals['n_interface_contacts'].append(float(nc))
            sample_vals['avg_interface_plddt'].append(aip)

    if not sample_conf_files and not sample_cif_files:
        print(f"WARNING: no confidence JSON or CIF for {pair_id}", file=sys.stderr)

    # =====================================================================
    # Compute mean and std, write TSV (same format as Boltz-2 metrics)
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
          f"iptm = {np.mean(sample_vals['iptm']):.4f} "
          f"+/- {np.std(sample_vals['iptm'], ddof=0):.4f}"
          if sample_vals['iptm'] else f"{pair_id}: no confidence data",
          file=sys.stderr)
    """
}
