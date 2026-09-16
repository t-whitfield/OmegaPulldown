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
    tuple val(pair_id), val(bait_id), val(bait_seq), val(candidate_id), val(candidate_seq), path(bait_msa), path(candidate_msa)
    val bait_ptms_path
    val candidate_ptms_path

    output:
    tuple val(pair_id), path("${pair_id}.json"), path("msa_out/*.a3m"), emit: json_with_msa

    script:
    def bait_ptms_arg = (bait_ptms_path != 'no_file') ? "--bait-ptms ${bait_ptms_path}" : ""
    def cand_ptms_arg = (candidate_ptms_path != 'no_file') ? "--candidate-ptms ${candidate_ptms_path}" : ""
    """
    make_af3_json.py \
        --bait-id "${bait_id}" \
        --bait-seq "${bait_seq}" \
        --bait-msa "${bait_msa}" \
        --candidate-id "${candidate_id}" \
        --candidate-seq "${candidate_seq}" \
        --candidate-msa "${candidate_msa}" \
        --seeds "${params.af3_seeds}" \
        ${bait_ptms_arg} \
        ${cand_ptms_arg} \
        --output "${pair_id}.json"

    # Collect AF3-ready A3M files into msa_out/.
    # make_af3_json.py creates *_af3.a3m files (cleaned + query-patched).
    # Using a new directory ensures Nextflow captures them as outputs.
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

    singularity run \
        --nv \
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
}
