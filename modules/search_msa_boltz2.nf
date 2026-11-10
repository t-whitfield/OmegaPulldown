/*
 * Boltz-2 local MSA search + taxonomy tagging (predictor=boltz2, local mode
 * only -- not used when params.use_msa_server=true).
 *
 * Split into three processes:
 *
 *   1. MSA_SEARCH_DB_BOLTZ2  -- UniRef30 search only (process_mmseqs_large,
 *      240GB). Always runs.
 *   2. MSA_SEARCH_ENV_BOLTZ2 -- env/metagenomic search only
 *      (process_mmseqs_bigmem, 800GB). Only invoked when params.use_env=true
 *      (see main.nf): this is the only piece of the chain that needs a
 *      typlically scarce bigmem-class node.
 *   3. MERGE_TAXONOMY_BOLTZ2 -- merges UniRef30 + env a3ms (or passes
 *      UniRef30 through untouched when env didn't run) and does the
 *      taxonomy-ID tagging + CSV write (process_medium, 128GB).
 */

/*
 * UniRef30 search only: the real `mmseqs_search_monomer` recipe.
 *
 * Batched across proteins, same rationale and mechanism as AF3-native's
 * MSA_SEARCH_DB (see search_msa_af3.nf): the expensive part is loading the
 * large UniRef30 index, not the search itself for one short query, so
 * batching amortizes that cost across all proteins in the batch.
 *
 * Emits qdb (the raw query DB: createdb output, several small files) and
 * prof_res (the profile built by the UniRef iterative search) in addition to
 * uniref.a3m: MSA_SEARCH_ENV_BOLTZ2 needs both as its own query/context if
 * use_env=true, and MERGE_TAXONOMY_BOLTZ2 needs qdb.lookup regardless (to
 * recover prot_id at unpackdb time).
 *
 * Label: process_mmseqs_large (240GB), same tier AF3-native's own
 * MSA_SEARCH_DB uses: NOT process_mmseqs_bigmem.
 */
process MSA_SEARCH_DB_BOLTZ2 {
    tag "batch_of_${batch.size()}"
    label 'process_mmseqs_large'

    input:
    val(batch)

    output:
    tuple val(batch), path("qdb*"), path("prof_res*"), path("uniref.a3m*"), emit: uniref_result

    script:
    def fasta_lines = batch.collect { prot_id, sequence -> ">${prot_id}\n${sequence}" }.join('\n')
    def db = params.boltz2_uniref30_db
    """
    eval "\$(${params.conda_init} shell.bash hook)" && conda activate ${params.taxonomy_env}

    # Write query FASTA -- one record per protein in this batch
    cat << 'EOF' > query.fasta
${fasta_lines}
EOF

    mmseqs createdb query.fasta qdb

    # Faithful reproduction of colabfold's mmseqs_search_monomer()
    # (filter=True branch, UniRef30 half).
    #
    # --db-load-mode 3 ("mmap+touch") on every subcommand, not colabfold's
    # own default of 2 ("mmap"): mode 2's lazy  mmap can stall for a very
    # long time, independent of available memory.
    mmseqs search qdb ${db} res tmp \
        --threads ${task.cpus} \
        --num-iterations ${params.mmseqs_iterations} \
        --db-load-mode 3 \
        -a -s ${params.boltz2_mmseqs_sensitivity} -e ${params.boltz2_mmseqs_evalue} --max-seqs 10000

    mmseqs expandaln qdb ${db}.idx res ${db}.idx res_exp \
        --db-load-mode 3 --threads ${task.cpus} \
        --expansion-mode 0 -e inf --expand-filter-clusters 1 --max-seq-id 0.95

    mmseqs mvdb tmp/latest/profile_1 prof_res
    mmseqs lndb qdb_h prof_res_h

    mmseqs align prof_res ${db}.idx res_exp res_exp_realign \
        --db-load-mode 3 -e 10 --max-accept 100000 --threads ${task.cpus} --alt-ali 10 -a

    mmseqs filterresult qdb ${db}.idx res_exp_realign res_exp_realign_filter \
        --db-load-mode 3 --qid 0 --qsc 0.8 --diff 0 --threads ${task.cpus} \
        --max-seq-id 1.0 --filter-min-enable 100

    mmseqs result2msa qdb ${db}.idx res_exp_realign_filter uniref.a3m \
        --msa-format-mode 6 --db-load-mode 3 --threads ${task.cpus} \
        --filter-msa 1 --filter-min-enable 1000 --diff 3000 \
        --qid 0.0,0.2,0.4,0.6,0.8,1.0 --qsc 0 --max-seq-id 0.95

    # "res*" matches only the search/expandaln/align/filterresult
    # intermediates (res, res_exp, res_exp_realign, res_exp_realign_filter),
    # NOT "prof_res*" (doesn't start with "res") or "qdb*"/"uniref.a3m*",
    # which are declared outputs and must survive this task.
    rm -rf tmp res*
    """

    stub:
    // No mmseqs invocation -- avoids real compute against production DBs
    // under -stub-run. Minimal stand-ins for all three declared outputs so
    // downstream processes (MSA_SEARCH_ENV_BOLTZ2, MERGE_TAXONOMY_BOLTZ2)
    // still see a valid shape under -stub-run.
    """
    eval "\$(${params.conda_init} shell.bash hook)" && conda activate ${params.taxonomy_env}
    cat << 'EOF' > query.fasta
${batch.collect { prot_id, sequence -> ">${prot_id}\n${sequence}" }.join('\n')}
EOF
    mmseqs createdb query.fasta qdb
    mmseqs createdb query.fasta prof_res
    mmseqs createdb query.fasta uniref.a3m
    """
}

/*
 * Env/metagenomic search only (use_env branch of colabfold's
 * mmseqs_search_monomer().
 *
 * Query is prof_res (the profile from MSA_SEARCH_DB_BOLTZ2's UniRef
 * iterative search), not a fresh qdb query: the whole point of running
 * this second, profile-seeded search rather than independently from
 * scratch, the same as colabfold's recipe.
 *
 * Label: process_mmseqs_bigmem (800GB).
 */
process MSA_SEARCH_ENV_BOLTZ2 {
    tag "batch_of_${batch.size()}"
    label 'process_mmseqs_bigmem'

    input:
    tuple val(batch), path(qdb_files), path(prof_res_files), path(uniref_a3m_files)

    output:
    tuple val(batch), path(qdb_files), path(uniref_a3m_files), path("env.a3m*"), emit: env_result

    script:
    def envdb = params.boltz2_envdb_db
    """
    eval "\$(${params.conda_init} shell.bash hook)" && conda activate ${params.taxonomy_env}

    mmseqs search prof_res ${envdb} res_env tmp3 \
        --threads ${task.cpus} \
        --num-iterations ${params.mmseqs_iterations} \
        --db-load-mode 3 \
        -a -s ${params.boltz2_mmseqs_sensitivity} -e ${params.boltz2_mmseqs_evalue} --max-seqs 10000

    mmseqs expandaln prof_res ${envdb}.idx res_env ${envdb}.idx res_env_exp \
        --db-load-mode 3 --threads ${task.cpus} \
        --expansion-mode 0 -e inf --expand-filter-clusters 1 --max-seq-id 0.95

    mmseqs align tmp3/latest/profile_1 ${envdb}.idx res_env_exp res_env_exp_realign \
        --db-load-mode 3 -e 10 --max-accept 100000 --threads ${task.cpus} --alt-ali 10 -a

    mmseqs filterresult qdb ${envdb}.idx res_env_exp_realign res_env_exp_realign_filter \
        --db-load-mode 3 --qid 0 --qsc 0.8 --diff 0 --threads ${task.cpus} \
        --max-seq-id 1.0 --filter-min-enable 100

    mmseqs result2msa qdb ${envdb}.idx res_env_exp_realign_filter env.a3m \
        --msa-format-mode 6 --db-load-mode 3 --threads ${task.cpus} \
        --filter-msa 1 --filter-min-enable 1000 --diff 3000 \
        --qid 0.0,0.2,0.4,0.6,0.8,1.0 --qsc 0 --max-seq-id 0.95

    # prof_res is only needed as this search's own query: not used by
    # MERGE_TAXONOMY_BOLTZ2, so it's fine to drop here rather than pass it
    # along further. "res_env*" matches only this process's own
    # search/expandaln/align/filterresult intermediates.
    rm -rf tmp3 res_env* prof_res*
    """

    stub:
    // Same minimal-stand-in approach as MSA_SEARCH_DB_BOLTZ2's own stub:
    // a fresh mmseqs DB under the expected output name, no real search.
    """
    eval "\$(${params.conda_init} shell.bash hook)" && conda activate ${params.taxonomy_env}
    cat << 'EOF' > query.fasta
${batch.collect { prot_id, sequence -> ">${prot_id}\n${sequence}" }.join('\n')}
EOF
    mmseqs createdb query.fasta env.a3m
    """
}

/*
 * Merge UniRef30 + env a3ms (when env ran) or pass UniRef30 through
 * untouched (when it didn't), then taxonomy-ID tag every row and write the
 * final Boltz-2 CSV directly (bin/a3m_to_boltz_csv.py), via a lookup against
 * Boltz's own real taxonomy.rdb, loaded into a local (per-batch, per-task)
 * redis-server: not a header-text regex (that information isn't present
 * in mmseqs2's own --msa-format-mode 6 output at all).
 *
 * Unlike AF3-native's MERGE_AF3_NATIVE_MSA (process_low, pure Python, no DB
 * access), this merge step is NOT trivial: it still does the full
 * taxonomy.rdb/redis load (~300s for 412M keys) and per-protein CSV
 * conversion, so it keeps the process_medium tier (128GB), not process_low.
 */
process MERGE_TAXONOMY_BOLTZ2 {
    tag "batch_of_${batch.size()}"
    label 'process_medium'

    publishDir "${params.outdir}/alignments_boltz2/csv", mode: 'copy'

    input:
    tuple val(batch), path(qdb_files), path(uniref_a3m_files), path(env_a3m_files)

    output:
    path("*.csv"), emit: csv_files

    script:
    """
    eval "\$(${params.conda_init} shell.bash hook)" && conda activate ${params.taxonomy_env}

    if [ -e env.a3m ]; then
        echo "=== use_env=true: merging UniRef30 + env a3ms ==="
        mmseqs mergedbs qdb final.a3m uniref.a3m env.a3m
    else
        echo "=== use_env=false: UniRef30 only, no merge needed ==="
        mmseqs mvdb uniref.a3m final.a3m
    fi

    # Split the packed result back into one file per protein (mergedbs/mvdb
    # output doesn't inherit qdb's .lookup automatically: same fix as
    # AF3-native's MSA_SEARCH_DB).
    cp qdb.lookup final.a3m.lookup
    mmseqs unpackdb final.a3m . --unpack-name-mode 1 --unpack-suffix ".a3m"

    rm -rf qdb* uniref.a3m* final.a3m* env.a3m* NO_ENV

    # Start a local (this task only, localhost, no cross-node reachability
    # needed), per-batch redis-server loaded with taxonomy.rdb.
    #
    # Cleanup via `trap ... EXIT`, not just a shutdown line at the end of the
    # script: a plain end-of-script shutdown only runs if execution reaches
    # it normally: if mmseqs2/python crash, the task hits this pipeline's
    # errorStrategy/maxRetries kill, or anything else exits the script early
    # (Nextflow wraps process scripts with `set -e` by default, so any
    # failing command exits immediately), redis-server would otherwise be
    # left running on the compute node indefinitely, holding ~16-27GB RAM and
    # port 7777. `trap ... EXIT` fires on every exit path, not just the happy
    # one.
    cp ${params.taxonomy_rdb} \$(basename ${params.taxonomy_rdb})
    redis-server --dir \$(pwd) --dbfilename \$(basename ${params.taxonomy_rdb}) --port 7777 &
    REDIS_PID=\$!
    # Capture the real exit status before cleanup and re-exit with it
    # explicitly, and keep `kill ... || true`: both fixes are essential,
    trap 'ec=\$?; redis-cli -p 7777 shutdown nosave >/dev/null 2>&1; kill \$REDIS_PID 2>/dev/null || true; exit \$ec' EXIT
    for i in \$(seq 1 120); do
        redis-cli -p 7777 ping 2>/dev/null | grep -q PONG && break
        sleep 5
    done

    for f in *.a3m; do
        prot_id="\${f%.a3m}"
        a3m_to_boltz_csv.py "\$f" --output "\${prot_id}.csv" \
            --max-seqs ${params.max_msa_seqs} --redis-port 7777
    done

    rm -f \$(basename ${params.taxonomy_rdb}) *.a3m
    """

    stub:
    // No mmseqs/redis invocation: avoids real compute against production
    // DBs and a 7.3GB taxonomy.rdb copy under -stub-run. One query-only CSV
    // (key=0 row only) per protein in the batch, a valid minimal MSA.
    def stub_lines = batch.collect { prot_id, sequence ->
        "printf 'sequence,key\\n${sequence},0\\n' > ${prot_id}.csv"
    }.join('\n')
    """
    ${stub_lines}
    """
}
