/*
 * AF3-native MSA search (predictor=alphafold3 only).
 *
 * Replicates AlphaFold 3's own genetic search using mmseqs2 as a validated
 * drop-in for jackhmmer/HHblits.
 *
 * Two independent outputs per unique protein, matching AF3-native's own
 * architecture:
 *   - Unpaired MSA: UniRef90 + MGnify + Small BFD, merged with
 *     exact-duplicate-sequence removal (MSA_SEARCH_DB + MERGE_AF3_NATIVE_MSA).
 *   - Paired MSA: UniProt only, kept completely separate, full original
 *     headers preserved so AlphaFold 3's own species-pairing code
 *     (alphafold3/data/msa_identifiers.py) can extract per-hit species IDs
 *     and do bait/candidate cross-chain matching itself at inference time
 *     (MSA_SEARCH_UNIPROT_PAIRED). No pairing logic is implemented here:
 *     AF3 does it internally given a correctly-headered pairedMsaPath.
 */

/*
 * db is a map: [name: <string>, path: <mmseqs2 db path>, max_seqs: <int>]
 * (see params.af3_native_dbs in nextflow.config). One invocation per
 * (batch-of-proteins, database): each database searched independently
 * against the raw queries, unlike ColabFold's cascaded profile search; no
 * profile is carried between databases, matching AF3-native's own
 * methodology.
 *
 * Batched across proteins: the expensive part of this process is loading
 * the large mmseqs2 index, not the search itself for one short
 * query: batching many proteins into one queryDB amortizes that load cost.
 * The multi-query  result is split back into per-protein files with
 * `mmseqs unpackdb`.
 */
process MSA_SEARCH_DB {
    tag "batch_of_${batch.size()}:${db.name}"
    // Static label, generously sized for all 3 DBs:
    // Nextflow's `label` directive doesn't support dynamic closures (unlike
    // memory/cpus) -- `label { db.label }` throws a cast exception. Small
    // BFD doesn't strictly need 240GB, but a single static over-provisioned
    // label is simpler and more robust than per-task label selection; the
    // only cost is some wasted headroom on BFD's own task, not a correctness
    // issue.
    label 'process_mmseqs_large'

    input:
    tuple val(batch), val(db)

    output:
    path("*.${db.name}.a3m"), emit: a3m_files

    script:
    // batch is a list of [prot_id, sequence] pairs (see main.nf's
    // split_ch.compute.collate(batch_size) -> combine(db_ch)).
    def fasta_lines = batch.collect { prot_id, sequence -> ">${prot_id}\n${sequence}" }.join('\n')
    def split_mem = "${task.memory.toGiga() - 10}G"
    """
    # Write query FASTA -- one record per protein in this batch
    cat << 'EOF' > query.fasta
${fasta_lines}
EOF

    # Create query database (one multi-record DB for the whole batch)
    mmseqs createdb query.fasta queryDB

    # Search against this single AF3-native database.
    # Flags match AlphaFast's own validated mmseqs2 recipe exactly: a
    # single-pass search (no --num-iterations; AF3-native genetic search is
    # not an iterative profile search), with -a for alignment backtraces
    # (required for result2msa to reconstruct the gapped alignment).
    mmseqs search queryDB ${db.path} resultDB tmp \
        --split-memory-limit ${split_mem} \
        --split ${db.split} \
        --split-mode 0 \
        -a \
        -s ${params.mmseqs_sensitivity} \
        -e ${params.af3_native_evalue} \
        --max-seqs ${db.max_seqs} \
        --threads ${task.cpus}

    # Convert to a3m MSA format (one packed multi-query result)
    mmseqs result2msa queryDB ${db.path} resultDB packed.a3m \
        --msa-format-mode 6

    # Split the packed result back into one file per protein.
    cp queryDB.lookup packed.a3m.lookup
    mmseqs unpackdb packed.a3m . --unpack-name-mode 1 --unpack-suffix ".${db.name}.a3m"

    # Clean up tmp
    rm -rf tmp queryDB* resultDB* packed.a3m packed.a3m.lookup
    """

    stub:
    // No mmseqs invocation: avoids real compute against production DBs
    // under -stub-run. One query-only a3m (no hits) per protein in the
    // batch, a valid minimal MSA: keeps -stub-run meaningful for the
    // batched shape, not just the single-protein case.
    def stub_lines = batch.collect { prot_id, sequence ->
        "echo '>${prot_id}' > ${prot_id}.${db.name}.a3m\necho '${sequence}' >> ${prot_id}.${db.name}.a3m"
    }.join('\n')
    """
    ${stub_lines}
    """
}

/*
 * Merge the three per-database a3m files (UniRef90, MGnify, Small BFD) for
 * one protein into a single unpaired MSA, with exact-duplicate-sequence
 * removal (see merge_af3_native_msa.py for the AF3-native dedup semantics).
 * UniProt is NOT one of the inputs (see MSA_SEARCH_UNIPROT_PAIRED).
 */
process MERGE_AF3_NATIVE_MSA {
    tag "${prot_id}"
    label 'process_low'

    publishDir "${params.outdir}/alignments_af3/a3m", mode: 'copy'

    input:
    tuple val(prot_id), path(a3m_files)

    output:
    tuple val(prot_id), path("${prot_id}.a3m"), emit: a3m

    script:
    """
    merge_af3_native_msa.py ${a3m_files} --output ${prot_id}.a3m
    """

    stub:
    // Pure-Python, local files only, no heavy/external dependency -- safe
    // and fast to run for real even under -stub-run.
    """
    merge_af3_native_msa.py ${a3m_files} --output ${prot_id}.a3m
    """
}

/*
 * UniProt search for the paired MSA, batched across proteins (see
 * MSA_SEARCH_DB above for the batching rationale).
 */
process MSA_SEARCH_UNIPROT_PAIRED {
    tag "batch_of_${batch.size()}:uniprot_paired"
    label 'process_medium'

    publishDir "${params.outdir}/alignments_af3/paired_a3m", mode: 'copy'

    input:
    val(batch)

    output:
    // Named "<prot_id>.paired.a3m", not "<prot_id>.a3m": these files get
    // staged into the same MAKE_AF3_JSON task directory as the unpaired
    // "<prot_id>.a3m" for the same protein: identical basenames would
    // collide when Nextflow stages both inputs together.
    path("*.paired.a3m"), emit: a3m_files

    script:
    def fasta_lines = batch.collect { prot_id, sequence -> ">${prot_id}\n${sequence}" }.join('\n')
    def split_mem = "${task.memory.toGiga() - 10}G"
    """
    # Write query FASTA -- one record per protein in this batch
    cat << 'EOF' > query.fasta
${fasta_lines}
EOF

    # Create query database (one multi-record DB for the whole batch)
    mmseqs createdb query.fasta queryDB

    mmseqs search queryDB ${params.af3_native_uniprot_db.path} resultDB tmp \
        --split-memory-limit ${split_mem} \
        --split ${params.af3_native_uniprot_db.split} \
        --split-mode 0 \
        -a \
        -s ${params.mmseqs_sensitivity} \
        -e ${params.af3_native_evalue} \
        --max-seqs ${params.af3_native_uniprot_db.max_seqs} \
        --threads ${task.cpus}

    # Aligned FASTA (mode 2), not A3M (mode 6) -- preserves full headers.
    # One packed multi-query result for the whole batch.
    mmseqs result2msa queryDB ${params.af3_native_uniprot_db.path} resultDB packed_aligned.fasta \
        --msa-format-mode 2

    # Split the packed result back into one aligned-FASTA file per protein.
    cp queryDB.lookup packed_aligned.fasta.lookup
    mmseqs unpackdb packed_aligned.fasta . --unpack-name-mode 1 --unpack-suffix ".aligned.fasta"

    # Per-protein a3m conversion: cheap pure-Python step, done in a loop
    # here rather than fanning back out to one task per protein.
    for f in *.aligned.fasta; do
        prot_id="\${f%.aligned.fasta}"
        aligned_fasta_to_a3m.py "\$f" --output "\${prot_id}.paired.a3m"
    done

    # Clean up tmp
    rm -rf tmp queryDB* resultDB* packed_aligned.fasta packed_aligned.fasta.lookup *.aligned.fasta
    """

    stub:
    // No mmseqs invocation: avoids real compute against production DBs
    // under -stub-run. One query-only a3m (no hits) per protein in the
    // batch, a valid minimal MSA.
    def stub_lines = batch.collect { prot_id, sequence ->
        "echo '>${prot_id}' > ${prot_id}.paired.a3m\necho '${sequence}' >> ${prot_id}.paired.a3m"
    }.join('\n')
    """
    ${stub_lines}
    """
}
