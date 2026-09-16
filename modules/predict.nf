/*
 * Run Boltz-2 structure prediction for each protein pair.
 *
 * Each YAML is predicted independently, enabling parallelism
 * across available GPUs via Nextflow's resource scheduler.
 *
 * Output: CIF structure, confidence JSON, PAE/PDE/pLDDT numpy arrays
 */

process BOLTZ_PREDICT {
    tag "${pair_id}"
    label 'process_gpu'

    beforeScript "eval \"\$(${params.conda_init} shell.bash hook)\" && conda activate ${params.boltz_env}"

    publishDir "${params.outdir}/structures_boltz2/${pair_id}", mode: 'copy',
        saveAs: { filename -> file(filename).name }

    input:
    tuple val(pair_id), path(yaml_file), path(msa_csvs)

    output:
    tuple val(pair_id), path("boltz_out/boltz_results_${pair_id}/predictions/${pair_id}/*"), emit: results

    script:
    def msa_flag      = params.use_msa_server ? "--use_msa_server" : ""
    def msa_url       = (params.use_msa_server && params.msa_server_url) \
                        ? "--msa_server_url ${params.msa_server_url}" : ""
    def potentials    = params.use_potentials ? "--use_potentials" : ""
    """
    boltz predict ${yaml_file} \
        --out_dir boltz_out \
        --recycling_steps ${params.recycling_steps} \
        --sampling_steps ${params.sampling_steps} \
        --diffusion_samples ${params.diffusion_samples} \
        --step_scale ${params.step_scale} \
        ${potentials} \
        ${msa_flag} \
        ${msa_url}
    """
}
