/*
 * Prepare Boltz-2 YAML input files for each bait x candidate pair.
 *
 * Takes paired protein sequences and pre-computed MSA CSV files,
 * generates a YAML file per pair in Boltz-2 input format.
 *
 * Optional PTM CSV files (protein_id,position,ccd) are passed through
 * to make_boltz_yaml.py, which filters rows by protein ID per pair.
 */

process GENERATE_PAIRS {
    label 'process_low'

    publishDir "${params.outdir}/pairs", mode: 'copy'

    input:
    path baits_fasta
    path candidates_fasta

    output:
    path "pairs.tsv", emit: pairs_tsv

    script:
    def self_flag = params.self_pairs ? "--self-pairs" : ""
    """
    generate_pairs.py \
        --baits ${baits_fasta} \
        --candidates ${candidates_fasta} \
        --output pairs.tsv \
        ${self_flag}
    """
}

process MAKE_BOLTZ_YAML {
    tag "${pair_id}"
    label 'process_low'

    publishDir "${params.outdir}/yamls", mode: 'copy'

    input:
    tuple val(pair_id), val(bait_id), val(bait_seq), val(candidate_id), val(candidate_seq), path(bait_msa), path(candidate_msa)
    val bait_ptms_path
    val candidate_ptms_path

    output:
    tuple val(pair_id), path("${pair_id}.yaml"), path(bait_msa), path(candidate_msa), emit: yaml_with_msa

    script:
    def bait_ptms_arg = (bait_ptms_path != 'no_file') ? "--bait-ptms ${bait_ptms_path}" : ""
    def cand_ptms_arg = (candidate_ptms_path != 'no_file') ? "--candidate-ptms ${candidate_ptms_path}" : ""
    """
    make_boltz_yaml.py \
        --bait-id "${bait_id}" \
        --bait-seq "${bait_seq}" \
        --bait-msa "${bait_msa}" \
        --candidate-id "${candidate_id}" \
        --candidate-seq "${candidate_seq}" \
        --candidate-msa "${candidate_msa}" \
        ${bait_ptms_arg} \
        ${cand_ptms_arg} \
        --output "${pair_id}.yaml"
    """
}

/*
 * Alternative: generate YAML without pre-computed MSA.
 * Boltz-2 will fetch MSAs from the server at prediction time.
 * Simpler but recomputes MSAs for shared proteins.
 */
process MAKE_BOLTZ_YAML_NO_MSA {
    tag "${pair_id}"
    label 'process_low'

    publishDir "${params.outdir}/yamls", mode: 'copy'

    input:
    tuple val(pair_id), val(bait_id), val(bait_seq), val(candidate_id), val(candidate_seq)
    val bait_ptms_path
    val candidate_ptms_path

    output:
    tuple val(pair_id), path("${pair_id}.yaml"), emit: yaml

    script:
    def bait_ptms_arg = (bait_ptms_path != 'no_file') ? "--bait-ptms ${bait_ptms_path}" : ""
    def cand_ptms_arg = (candidate_ptms_path != 'no_file') ? "--candidate-ptms ${candidate_ptms_path}" : ""
    """
    make_boltz_yaml.py \
        --bait-id "${bait_id}" \
        --bait-seq "${bait_seq}" \
        --candidate-id "${candidate_id}" \
        --candidate-seq "${candidate_seq}" \
        ${bait_ptms_arg} \
        ${cand_ptms_arg} \
        --output "${pair_id}.yaml"
    """
}
