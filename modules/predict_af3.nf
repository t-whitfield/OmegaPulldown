/*
 * AlphaFold 3 structure prediction processes.
 *
 * MAKE_AF3_JSON: Generate AF3-format JSON input with pre-computed MSAs.
 * AF3_PREDICT:   Run AlphaFold 3 via Singularity container.
 *
 * MSA files (A3M) are staged alongside the JSON so that AF3 can
 * resolve relative unpairedMsaPath references.  The --norun_data_pipeline
 * flag skips jackhmmer/nhmmer searches, using the provided MSAs directly.
 */

process MAKE_AF3_JSON {
    tag "${pair_id}"
    label 'process_low'

    publishDir "${params.outdir}/jsons", mode: 'copy'

    input:
    tuple val(pair_id), val(bait_id), val(bait_seq), val(candidate_id), val(candidate_seq), path(bait_msa), path(candidate_msa), path(bait_paired_msa), path(candidate_paired_msa)
    val ptms_path

    output:
    tuple val(pair_id), path("${pair_id}.json"), path("msa_out/*.a3m"), emit: json_with_msa

    script:
    // One ID-keyed PTM table (or a comma-separated list of them); the
    // script resolves each chain's PTMs by protein ID, never by bait or
    // candidate role. See README.pair_input.md.
    def ptms_arg = (ptms_path != 'no_file') ? "--ptms ${ptms_path}" : ""
    """
    make_af3_json.py \
        --bait-id "${bait_id}" \
        --bait-seq "${bait_seq}" \
        --bait-msa "${bait_msa}" \
        --bait-paired-msa "${bait_paired_msa}" \
        --candidate-id "${candidate_id}" \
        --candidate-seq "${candidate_seq}" \
        --candidate-msa "${candidate_msa}" \
        --candidate-paired-msa "${candidate_paired_msa}" \
        --seeds "${params.af3_seeds}" \
        ${ptms_arg} \
        --output "${pair_id}.json"

    # Collect AF3-ready A3M files into msa_out/.
    # make_af3_json.py creates *_af3.a3m files (cleaned + query-patched).
    # Using a new directory ensures Nextflow captures them as outputs.
    mkdir -p msa_out
    cp -L *_af3.a3m msa_out/
    """

    stub:
    // Pure-Python, local files only: safe and fast to run for real.
    // One ID-keyed PTM table (or a comma-separated list of them); the
    // script resolves each chain's PTMs by protein ID, never by bait or
    // candidate role. See README.pair_input.md.
    def ptms_arg = (ptms_path != 'no_file') ? "--ptms ${ptms_path}" : ""
    """
    make_af3_json.py \
        --bait-id "${bait_id}" \
        --bait-seq "${bait_seq}" \
        --bait-msa "${bait_msa}" \
        --bait-paired-msa "${bait_paired_msa}" \
        --candidate-id "${candidate_id}" \
        --candidate-seq "${candidate_seq}" \
        --candidate-msa "${candidate_msa}" \
        --candidate-paired-msa "${candidate_paired_msa}" \
        --seeds "${params.af3_seeds}" \
        ${ptms_arg} \
        --output "${pair_id}.json"

    mkdir -p msa_out
    cp -L *_af3.a3m msa_out/
    """
}

process AF3_PREDICT {
    tag "${pair_id}"
    label 'process_gpu'

    publishDir "${params.outdir}/structures_af3/${pair_id}", mode: 'copy',
        saveAs: { filename -> file(filename).name }

    input:
    tuple val(pair_id), path(json_file), path(msa_files)

    output:
    tuple val(pair_id), path("af3_out/${pair_id}/**"), emit: results

    script:
    """
    # Nextflow stages input files as symlinks, but Singularity bind mounts
    # cannot follow symlinks to paths outside the container.  Copy the
    # JSON and A3M files into a local staging directory.
    mkdir -p staging
    cp -L ${json_file} staging/
    for f in ${msa_files}; do
        cp -L \$f staging/
    done

    # Unified memory plus an enlarged XLA budget, so a complex whose working
    # set exceeds physical VRAM spills to host RAM instead of dying with
    # RESOURCE_EXHAUSTED.
    VRAM_MIB=\$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits | head -1)
    MEM_FRACTION=\$(awk -v a=${task.memory.toMega()} -v v=\$VRAM_MIB 'BEGIN{f=0.9*a/v; if(f<1)f=1; if(f>4)f=4; printf "%.2f", f}')
    echo "AF3 memory budget: VRAM \${VRAM_MIB} MiB, alloc ${task.memory.toMega()} MiB, XLA_CLIENT_MEM_FRACTION=\${MEM_FRACTION}"

    singularity run \
        --nv \
        --env TF_FORCE_UNIFIED_MEMORY=1,XLA_CLIENT_MEM_FRACTION=\${MEM_FRACTION} \
        -B \${PWD}/staging:/af3_work \
        -B ${params.af3_db_dir}:/af3_db \
        -B ${params.af3_model_dir}:/af3_models \
        ${params.af3_sif} \
        python /app/alphafold/run_alphafold.py \
            --json_path=/af3_work/${json_file} \
            --model_dir=/af3_models \
            --db_dir=/af3_db \
            --output_dir=/af3_work/af3_out \
            --norun_data_pipeline \
            --num_diffusion_samples=${params.af3_samples}

    # Move output back to work dir so Nextflow can capture it
    mv staging/af3_out af3_out
    """

    stub:
    // No singularity/GPU invocation: this needs real model weights and a
    // GPU, neither is appropriate for a stub-run on the head node. Minimal
    // placeholder output matching what PARSE_METRICS_AF3 (rank_af3.nf)
    // globs for: *_sample-*_summary_confidences.json, *_sample-*_model.cif.
    """
    mkdir -p af3_out/${pair_id}
    echo '{"ranking_score": 0.5, "ptm": 0.5, "iptm": 0.5}' > af3_out/${pair_id}/${pair_id}_seed-42_sample-0_summary_confidences.json
    touch af3_out/${pair_id}/${pair_id}_seed-42_sample-0_model.cif
    """
}
