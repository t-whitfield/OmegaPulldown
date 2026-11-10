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
    import json, glob, math, re, sys, os, shutil
    import numpy as np

    pair_id = "${pair_id}"

    # pDockQ comes from bin/compute_pdockq.py.
    _p = shutil.which("compute_pdockq.py")
    if _p is None:
        sys.exit("compute_pdockq.py not found on PATH -- is bin/ present?")
    sys.path.insert(0, os.path.dirname(_p))
    from compute_pdockq import compute_pdockq

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

    stub:
    // Minimal placeholder, not a duplicate of the embedded Python parser
    // above. This just satisfies the declared output for wiring checks.
    """
    printf "pair_id\\tn_samples\\tconfidence_score_mean\\tconfidence_score_std\\n" > ${pair_id}.metrics.tsv
    printf "${pair_id}\\t1\\t0.5\\tNA\\n" >> ${pair_id}.metrics.tsv
    """
}
