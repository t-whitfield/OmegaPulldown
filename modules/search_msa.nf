/*
 * MSA generation using MMseqs2 against a local sequence database.
 *
 * Computes one MSA per unique protein sequence (not per pair),
 * so shared baits/targets are searched only once. Outputs a3m files
 * which are subsequently converted to Boltz-2 CSV format.
 *
 * Requires: mmseqs2 in PATH, local sequence database
 *   (e.g., colabfold_envdb, uniref30)
 */

process MMSEQS_SEARCH {
    tag "${prot_id}"
    label 'process_medium'

    publishDir "${params.outdir}/alignments/a3m", mode: 'copy'

    input:
    tuple val(prot_id), val(sequence)

    output:
    tuple val(prot_id), path("${prot_id}.a3m"), emit: a3m

    script:
    """
    # Write query FASTA
    echo ">${prot_id}" > query.fasta
    echo "${sequence}" >> query.fasta

    # Create query database
    mmseqs createdb query.fasta queryDB

    # Search against target database
    # Note, the "split-memory-limit" setting is determined by
    # cluster capacity and utilization.
    mmseqs search queryDB ${params.mmseqs_db} resultDB tmp \
        --split-memory-limit 120G \
	--num-iterations ${params.mmseqs_iterations} \
        -s ${params.mmseqs_sensitivity} \
        -e ${params.mmseqs_evalue} \
        --threads ${task.cpus}

    # Convert to a3m MSA format
    mmseqs result2msa queryDB ${params.mmseqs_db} resultDB ${prot_id}.a3m \
        --msa-format-mode 6

    # Clean up tmp
    rm -rf tmp queryDB* resultDB*
    """
}

process A3M_TO_BOLTZ_CSV {
    tag "${prot_id}"
    label 'process_low'

    publishDir "${params.outdir}/alignments/csv", mode: 'copy'

    input:
    tuple val(prot_id), path(a3m_file)

    output:
    tuple val(prot_id), path("${prot_id}.csv"), emit: csv

    script:
    """
    a3m_to_boltz_csv.py ${a3m_file} --output ${prot_id}.csv \
        --max-seqs ${params.max_msa_seqs}
    """
}
