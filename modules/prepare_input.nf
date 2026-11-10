/*
 * Pair enumeration (predictor-agnostic) and Boltz-2 YAML assembly.
 *
 * GENERATE_PAIRS builds the pair list that drives the whole pipeline, in
 * either of two modes: cross-product (baits x candidates) or pair-list
 * (only the pairs named in params.pairs_tsv). Both emit the same TSV, so
 * every downstream stage is mode-agnostic; see bin/generate_pairs.py and
 * README.pair_input.md. Pairs are collapsed on tuple(sorted((a,b))) since
 * A_with_B and B_with_A are the same complex: a no-op when the bait and
 * candidate sets are disjoint, and a halving when they overlap.
 *
 * MAKE_BOLTZ_YAML / MAKE_BOLTZ_YAML_NO_MSA then assemble one Boltz-2 input
 * YAML per pair, with or without pre-computed MSA CSVs. (The AF3 equivalent,
 * MAKE_AF3_JSON, lives in predict_af3.nf.)
 *
 * PTMs come from a single id-keyed table (params.ptms, conventionally
 * ptms.tsv; tab or comma). make_boltz_yaml.py resolves each chain's
 * modifications by protein ID, never by whether that protein is the bait or
 * the candidate of this particular pair.
 */

process GENERATE_PAIRS {
    label 'process_low'

    publishDir "${params.outdir}/pairs", mode: 'copy'

    input:
    path baits_fasta, stageAs: 'pairgen_input_a'
    path candidates_fasta, stageAs: 'pairgen_input_b'

    output:
    path "pairs.tsv", emit: pairs_tsv

    script:
    // Two input modes, same output contract (see bin/generate_pairs.py).
    // Cross-product: baits x candidates, the classic pulldown design.
    // Pair-list: only the pairs named in params.pairs_tsv, for adding
    // structure-based ranking to an experimental edge list. main.nf decides
    // which files are staged here; this process only picks the matching
    // command line. --self-pairs applies to the cross product only: an
    // explicitly listed self-pair is always kept, since the list is the
    // user's stated intent.
    def self_flag = params.self_pairs ? "--self-pairs" : ""
    if (params.pairs_tsv)
        """
        generate_pairs.py \
            --pairs-tsv ${baits_fasta} \
            --sequences ${candidates_fasta} \
            --output pairs.tsv
        """
    else
        """
        generate_pairs.py \
            --baits ${baits_fasta} \
            --candidates ${candidates_fasta} \
            --output pairs.tsv \
            ${self_flag}
        """

    stub:
    // Pure-Python, local files only: safe and fast to run for real; also
    // gives downstream stub processes a genuine pairs.tsv to fan out over.
    def self_flag = params.self_pairs ? "--self-pairs" : ""
    if (params.pairs_tsv)
        """
        generate_pairs.py \
            --pairs-tsv ${baits_fasta} \
            --sequences ${candidates_fasta} \
            --output pairs.tsv
        """
    else
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
    val ptms_path

    output:
    tuple val(pair_id), path("${pair_id}.yaml"), path(bait_msa), path(candidate_msa), emit: yaml_with_msa

    script:
    // One ID-keyed PTM table (or a comma-separated list of them); the
    // script resolves each chain's PTMs by protein ID, never by bait or
    // candidate role.
    def ptms_arg = (ptms_path != 'no_file') ? "--ptms ${ptms_path}" : ""
    """
    make_boltz_yaml.py \
        --bait-id "${bait_id}" \
        --bait-seq "${bait_seq}" \
        --bait-msa "${bait_msa}" \
        --candidate-id "${candidate_id}" \
        --candidate-seq "${candidate_seq}" \
        --candidate-msa "${candidate_msa}" \
        ${ptms_arg} \
        --output "${pair_id}.yaml"
    """

    stub:
    // Pure-Python, local files only: safe and fast to run for real.
    // One ID-keyed PTM table (or a comma-separated list of them); the
    // script resolves each chain's PTMs by protein ID, never by bait or
    // candidate role.
    def ptms_arg = (ptms_path != 'no_file') ? "--ptms ${ptms_path}" : ""
    """
    make_boltz_yaml.py \
        --bait-id "${bait_id}" \
        --bait-seq "${bait_seq}" \
        --bait-msa "${bait_msa}" \
        --candidate-id "${candidate_id}" \
        --candidate-seq "${candidate_seq}" \
        --candidate-msa "${candidate_msa}" \
        ${ptms_arg} \
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
    val ptms_path

    output:
    tuple val(pair_id), path("${pair_id}.yaml"), emit: yaml

    script:
    // One ID-keyed PTM table (or a comma-separated list of them); the
    // script resolves each chain's PTMs by protein ID, never by bait or
    // candidate role.
    def ptms_arg = (ptms_path != 'no_file') ? "--ptms ${ptms_path}" : ""
    """
    make_boltz_yaml.py \
        --bait-id "${bait_id}" \
        --bait-seq "${bait_seq}" \
        --candidate-id "${candidate_id}" \
        --candidate-seq "${candidate_seq}" \
        ${ptms_arg} \
        --output "${pair_id}.yaml"
    """

    stub:
    // Pure-Python, local files only: safe and fast to run for real.
    // One ID-keyed PTM table (or a comma-separated list of them); the
    // script resolves each chain's PTMs by protein ID, never by bait or
    // candidate role.
    def ptms_arg = (ptms_path != 'no_file') ? "--ptms ${ptms_path}" : ""
    """
    make_boltz_yaml.py \
        --bait-id "${bait_id}" \
        --bait-seq "${bait_seq}" \
        --candidate-id "${candidate_id}" \
        --candidate-seq "${candidate_seq}" \
        ${ptms_arg} \
        --output "${pair_id}.yaml"
    """
}
