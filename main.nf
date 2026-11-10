/*
 * =========================================================================
 *  OmegaPulldown: Protein-protein interaction screen using Boltz-2
 *                 or AlphaFold 3
 * =========================================================================
 *
 *  Input:  either a bait FASTA + candidate FASTA (cross-product mode), or a
 *          list of pairs + one sequence FASTA (pair-list mode, --pairs_tsv).
 *  Output: ranked table of predicted complex structures with confidence
 *          metrics (ipTM, pDockQ, pLDDT, pTM), plus per-pair structure
 *          subdirectories under results/structures_boltz2/ or
 *          results/structures_af3/ depending on the predictor.
 *
 *  Pipeline stages:
 *    1. GENERATE_PAIRS - build the pair list. Cross-product mode enumerates
 *                        bait x candidate; pair-list mode takes only the
 *                        pairs named in --pairs_tsv. Either way pairs are
 *                        collapsed on tuple(sorted((a,b))), since A_with_B
 *                        and B_with_A are the same complex (a no-op when the
 *                        bait and candidate sets are disjoint).
 *
 *    2a. (Boltz-2 local MSA, use_msa_server=false)
 *        MSA_SEARCH_DB_BOLTZ2, UniRef30 search per unique protein (CPU);
 *        then (only if params.use_env=true) MSA_SEARCH_ENV_BOLTZ2,
 *        env/metagenomic search; then MERGE_TAXONOMY_BOLTZ2, merge +
 *        taxonomy-ID tagging into the final per-protein CSVs (see
 *        modules/search_msa_boltz2.nf). No separate pairing/join step:
 *        Boltz-2's own featurizer pairs rows sharing a taxonomy ID across
 *        whichever two chains' CSVs are loaded for a prediction.
 *
 *    2b. (AF3-native MSA, predictor=alphafold3)
 *        MSA_SEARCH_DB, one task per unpaired database (UniRef90, MGnify,
 *        Small BFD); then MERGE_AF3_NATIVE_MSA, merge with exact-duplicate
 *        removal into one unpaired a3m per protein. Separately,
 *        MSA_SEARCH_UNIPROT_PAIRED builds the PAIRED a3m from UniProt, which
 *        AF3 keeps distinct from the unpaired MSA and uses for cross-chain
 *        species matching (see modules/search_msa_af3.nf).
 *
 *    3. MAKE_BOLTZ_YAML / MAKE_AF3_JSON - assemble input per pair
 *    4. BOLTZ_PREDICT / AF3_PREDICT     - structure prediction (GPU)
 *    5. PARSE_METRICS / PARSE_METRICS_AF3 - confidence + pDockQ per pair
 *                        (pDockQ via bin/compute_pdockq.py, shared by both)
 *    6. RANK_PREDICTIONS - aggregate and rank all pairs
 *
 *  Predictor selection:
 *    --predictor boltz2       (default)
 *    --predictor alphafold3
 *
 *  MSA caching -- NOTE the two predictors cache different files:
 *    predictor=boltz2 caches per-protein CSVs at
 *      params.alignments_dir_boltz2/csv/<prot_id>.csv
 *    There is no a3m-level cache for Boltz-2: taxonomy tagging happens in
 *    the same task as the search, so no intermediate a3m is worth keeping.
 *    predictor=alphafold3 caches per-protein a3m files at
 *      params.alignments_dir_af3/a3m/<prot_id>.a3m
 *      params.alignments_dir_af3/paired_a3m/<prot_id>.paired.a3m
 *    A protein whose footprint already exists skips its search entirely.
 *    Boltz-2 and AF3 use separate alignment directories: they were trained
 *    on different MSA lineages, so a shared directory is not a neutral cache.
 *
 *  When params.use_msa_server = true, stage 2a is skipped and Boltz-2
 *  fetches MSAs from the ColabFold API at prediction time.
 *  (Only valid with predictor = boltz2.)
 *
 *  Post-translational modifications (optional):
 *    Place ptms.tsv in the same directory as the sequence FASTA, or specify
 *    it explicitly via --ptms. Modifications are resolved by protein ID, so
 *    one table covers both chains of every pair.  Format (tab or comma):
 *      protein_id	position	ccd
 *      BAIT1_GFP	18	SEP
 *
 *  Usage:
 *    # Boltz-2 (default), API-based MSA (use_msa_server=true is the default --
 *    # recommended unless running a large all-vs-all screen, see
 *    # boltz2_uniref30_db/taxonomy_rdb in nextflow.config for the local-MSA
 *    # alternative)
 *    nextflow run main.nf \
 *      --baits_fasta inputs/sequences/baits.fasta \
 *      --candidates_fasta inputs/sequences/candidates.fasta
 *
 *    # Boltz-2, local MSA mode (no public API dependency and recommended for
 *    # large all-vs-all screens; see modules/search_msa_boltz2.nf). UniRef30
 *    # search only by default (params.use_env=false, see nextflow.config
 *    # for why); add --use_env true --boltz2_envdb_db /path/to/envdb to
 *    # also run the env/metagenomic search, at the cost of a scarce
 *    # bigmem-class node allocation.
 *    nextflow run main.nf \
 *      --baits_fasta inputs/sequences/baits.fasta \
 *      --candidates_fasta inputs/sequences/candidates.fasta \
 *      --use_msa_server false \
 *      --boltz2_uniref30_db /path/to/uniref30_2302_db \
 *      --taxonomy_rdb /path/to/taxonomy.rdb
 *
 *    # AlphaFold 3 (use_msa_server is not supported here: AF3 always uses
 *    # the locally built AF3-native MSAs of stage 2b)
 *    nextflow run main.nf \
 *      --predictor alphafold3 \
 *      --use_msa_server false \
 *      --baits_fasta inputs/sequences/baits.fasta \
 *      --candidates_fasta inputs/sequences/candidates.fasta
 *
 *    # Pair-list mode: predict ONLY the pairs named in a file, e.g. to add a
 *    # structure-based ranking to PPIs reported by proteomics. The first two
 *    # columns are protein IDs (header auto-detected); --sequences_fasta must
 *    # hold every protein the list names. The --baits_fasta and
 *    # --candidates_fasta flags are ignored in this mode.
 *    nextflow run main.nf \
 *      --pairs_tsv ms_hits.tsv \
 *      --sequences_fasta all_proteins.fasta
 *
 *    # Reuse pre-computed alignments (skips the MSA stage for any protein
 *    # whose cache is already present, see "MSA caching" above)
 *    nextflow run main.nf \
 *      --predictor alphafold3 \
 *      --use_msa_server false \
 *      --alignments_dir_af3 results/alignments_af3
 * =========================================================================
 */

nextflow.enable.dsl = 2

// --- Module imports ---
include { MSA_SEARCH_DB_BOLTZ2; MSA_SEARCH_ENV_BOLTZ2; MERGE_TAXONOMY_BOLTZ2 }                 from './modules/search_msa_boltz2'
include { MSA_SEARCH_DB; MERGE_AF3_NATIVE_MSA; MSA_SEARCH_UNIPROT_PAIRED }                     from './modules/search_msa_af3'
include { GENERATE_PAIRS; MAKE_BOLTZ_YAML; MAKE_BOLTZ_YAML_NO_MSA }   from './modules/prepare_input'
include { BOLTZ_PREDICT }                                             from './modules/predict_boltz2'
include { PARSE_METRICS; RANK_PREDICTIONS }                           from './modules/rank'
include { MAKE_AF3_JSON; AF3_PREDICT }                                from './modules/predict_af3'
include { PARSE_METRICS_AF3 }                                         from './modules/rank_af3'


/*
 * Parse a multi-sequence FASTA file into a list of [id, sequence] pairs.
 * Used at workflow scope to identify unique proteins for MSA computation.
 */
def parseFasta(fastaPath) {
    def seqs = []
    def f = file(fastaPath)
    def header = null
    def seqLines = []
    f.readLines().each { line ->
        if (line.startsWith('>')) {
            if (header != null) {
                seqs << [header, seqLines.join('')]
            }
            header = line.substring(1).split(/\s+/)[0]
            seqLines = []
        } else if (line.trim()) {
            seqLines << line.trim()
        }
    }
    if (header != null) {
        seqs << [header, seqLines.join('')]
    }
    return seqs
}


/*
 * Resolve the PTM table for this run.
 *
 * ONE ID-keyed table, conventionally ptms.tsv:
 *   protein_id<TAB>position<TAB>ccd
 * Tab- or comma-separated. Resolution by protein ID alone -- never
 * by whether a protein is the bait or the candidate of a given pair.
 *
 * Priority:
 *   1. Explicit --ptms
 *   2. Convention: ptms.tsv, then ptms.csv, beside the FASTA that defines
 *      this run's proteins (sequences_fasta in pair-list mode, otherwise
 *      baits_fasta)
 *   3. No PTMs ('no_file' sentinel)
 */
def resolvePtms(primaryFasta) {
    // 1. Explicit
    if (params.ptms != null) {
        def f = file(params.ptms)
        if (!f.exists()) {
            error "ptms file not found: ${params.ptms}"
        }
        log.info "PTM table (explicit): ${f}"
        return f.toString()
    }

    // 2. Convention
    if (primaryFasta != null) {
        def dir = file(primaryFasta).parent
        for (name in ['ptms.tsv', 'ptms.csv']) {
            def p = dir.resolve(name)
            if (p.exists()) {
                log.info "PTM table (auto-detected): ${p}"
                return p.toString()
            }
        }
    }

    // 3. None
    return 'no_file'
}

workflow {

    // ---------------------------------------------------------------
    // Validate predictor choice
    // ---------------------------------------------------------------
    def predictor = params.predictor.toLowerCase()
    if (!(predictor in ['boltz2', 'alphafold3'])) {
        error "Invalid predictor '${params.predictor}'. Choose 'boltz2' or 'alphafold3'."
    }
    if (predictor == 'alphafold3' && params.use_msa_server) {
        error "use_msa_server is not supported with predictor=alphafold3 " +
              "(it fetches Boltz-2-lineage MSAs from the ColabFold API; AF3 " +
              "needs its own AF3-native MSAs). Pass --use_msa_server false. " +
              "Note -profile test sets it true, so that profile needs the " +
              "override explicitly."
    }

    // Backend-specific alignments directory: Boltz-2 and AF3 were trained
    // on different MSA lineages, so they don't share a cache.
    def alignments_dir = (predictor == 'alphafold3') ? params.alignments_dir_af3 : params.alignments_dir_boltz2

    log.info """
    =========================================================================
    OmegaPulldown
    =========================================================================
    Predictor    : ${predictor}
    Baits        : ${params.baits_fasta}
    Candidates   : ${params.candidates_fasta}
    Alignments   : ${alignments_dir}
    MSA server   : ${params.use_msa_server}
    Output       : ${params.outdir}
    =========================================================================
    """.stripIndent()

    // ---------------------------------------------------------------
    // Resolve optional PTM CSV files
    // ---------------------------------------------------------------
    // The "primary" FASTA is whichever file defines this run's protein set:
    // sequences_fasta in pair-list mode, baits_fasta otherwise.
    ptms_path = resolvePtms(params.pairs_tsv ? params.sequences_fasta : params.baits_fasta)

    // ---------------------------------------------------------------
    // Stage 1: Build the pair list
    //
    // Cross-product mode (default): every bait x candidate pair.
    // Pair-list mode (--pairs_tsv): only the pairs named in that file,
    // with sequences drawn from --sequences_fasta. Intended for adding a
    // structure-based ranking to an experimental edge list (e.g. PPIs
    // reported by proteomics) rather than enumerating a cross product.
    //
    // Both modes feed GENERATE_PAIRS the same two path inputs, so the
    // process signature and every downstream stage are unchanged.
    // ---------------------------------------------------------------
    if (params.pairs_tsv) {
        def pairs_file = file(params.pairs_tsv)
        if (!pairs_file.exists()) {
            error "pairs_tsv not found: ${params.pairs_tsv}"
        }
        if (!params.sequences_fasta) {
            error "--pairs_tsv requires --sequences_fasta (a FASTA holding every protein the pair list names)."
        }
        def seqs_file = file(params.sequences_fasta)
        if (!seqs_file.exists()) {
            error "sequences_fasta not found: ${params.sequences_fasta}"
        }
        log.info "Pair-list mode: ${pairs_file} against ${seqs_file}"
        baits_ch      = Channel.fromPath(params.pairs_tsv)
        candidates_ch = Channel.fromPath(params.sequences_fasta)
    }
    else {
        baits_ch      = Channel.fromPath(params.baits_fasta)
        candidates_ch = Channel.fromPath(params.candidates_fasta)
    }

    GENERATE_PAIRS(baits_ch, candidates_ch)

    // Parse TSV -> channel of (pair_id, bait_id, bait_seq, candidate_id, candidate_seq)
    pairs_ch = GENERATE_PAIRS.out.pairs_tsv
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(
            row.pair_id,
            row.bait_id,
            row.bait_seq,
            row.candidate_id,
            row.candidate_seq
        )}

    if (params.use_msa_server && predictor == 'boltz2') {

        // =============================================================
        // Server mode: skip MSA pre-computation (Boltz-2 only)
        // Boltz-2 will use --use_msa_server to fetch MSAs on the fly
        // =============================================================

        predict_ch = MAKE_BOLTZ_YAML_NO_MSA(
            pairs_ch,
            ptms_path
        ).yaml
            .map { pair_id, yaml_file -> tuple(pair_id, yaml_file, []) }

    } else {

        // =============================================================
        // Local MSA mode: pre-compute MSAs per unique protein
        // Avoids redundant searches when proteins appear in many pairs.
        // Skips proteins that already have A3M files in alignments_dir.
        // =============================================================

        // Collect all unique protein sequences across both files
        bait_seqs      = parseFasta(params.baits_fasta)
        candidate_seqs = parseFasta(params.candidates_fasta)
        all_seqs       = (bait_seqs + candidate_seqs).unique { it[0] }

        unique_proteins_ch = Channel.from(all_seqs)
            .map { entry -> tuple(entry[0], entry[1]) }

        // Stage 2: MSA search for proteins without cached MSAs (CPU-bound).
        // predictor=alphafold3 uses the AF3-native 4-DB search,
        // caching per-protein a3m files; predictor=boltz2 uses
        // the ColabFold-lineage single-DB search + taxonomy tagging,
        // caching per-protein CSV files directly.
        if (predictor == 'alphafold3') {
            // --- MSA caching: split into cached vs. need-to-compute ---
            unique_proteins_ch.branch {
                cached:  file("${alignments_dir}/a3m/${it[0]}.a3m").exists()
                compute: true
            }.set { split_ch }

            // Use cached A3M files directly
            cached_a3m = split_ch.cached.map { prot_id, seq ->
                tuple(prot_id, file("${alignments_dir}/a3m/${prot_id}.a3m"))
            }

            // Batch proteins needing search into fixed-size groups before
            // invoking MSA_SEARCH_DB.
            db_ch = Channel.from(params.af3_native_dbs)
            batched_compute_ch = split_ch.compute.collate(params.af3_native_batch_size).map { batch -> [batch] }
            search_input_ch = batched_compute_ch.combine(db_ch)

            MSA_SEARCH_DB(search_input_ch)
            per_protein_a3m_ch = MSA_SEARCH_DB.out.a3m_files
                .flatten()
                .map { f -> tuple(f.name.tokenize('.')[0], f) }
                .groupTuple()

            MERGE_AF3_NATIVE_MSA(per_protein_a3m_ch)
            computed_a3m = MERGE_AF3_NATIVE_MSA.out.a3m

            // Merge cached + newly computed A3M streams
            all_a3m = cached_a3m.mix(computed_a3m)

            // Paired MSA (UniProt only): independent cache check, since this
            // is a separate type from the unpaired a3m above.
            unique_proteins_ch.branch {
                cached:  file("${alignments_dir}/paired_a3m/${it[0]}.paired.a3m").exists()
                compute: true
            }.set { paired_split_ch }

            cached_paired_msa = paired_split_ch.cached.map { prot_id, seq ->
                tuple(prot_id, file("${alignments_dir}/paired_a3m/${prot_id}.paired.a3m"))
            }

            // Same batching as the unpaired search above.
            batched_paired_compute_ch = paired_split_ch.compute.collate(params.af3_native_batch_size)
            MSA_SEARCH_UNIPROT_PAIRED(batched_paired_compute_ch)

            computed_paired_msa = MSA_SEARCH_UNIPROT_PAIRED.out.a3m_files
                .flatten()
                .map { f -> tuple(f.name.tokenize('.')[0], f) }

            all_paired_msa = cached_paired_msa.mix(computed_paired_msa)
        } else {
            // predictor == 'boltz2', local MSA mode (use_msa_server=false)

            // --- CSV caching: split into cached vs. need-to-compute ---
            // (CSV, not a3m, see the comment above this if/else for why
            // this predictor caches at the CSV level directly.)
            unique_proteins_ch.branch {
                cached:  file("${alignments_dir}/csv/${it[0]}.csv").exists()
                compute: true
            }.set { csv_split_ch }

            cached_csv = csv_split_ch.cached.map { prot_id, seq ->
                tuple(prot_id, file("${alignments_dir}/csv/${prot_id}.csv"))
            }

            batched_compute_ch = csv_split_ch.compute.collate(params.boltz2_batch_size)
            MSA_SEARCH_DB_BOLTZ2(batched_compute_ch)

            if (params.use_env) {
                MSA_SEARCH_ENV_BOLTZ2(MSA_SEARCH_DB_BOLTZ2.out.uniref_result)
                merge_input_ch = MSA_SEARCH_ENV_BOLTZ2.out.env_result
            } else {
                no_env_file = file("${projectDir}/assets/NO_ENV")
                merge_input_ch = MSA_SEARCH_DB_BOLTZ2.out.uniref_result
                    .map { batch, qdb, prof_res, uniref_a3m -> tuple(batch, qdb, uniref_a3m, no_env_file) }
            }
            MERGE_TAXONOMY_BOLTZ2(merge_input_ch)

            computed_csv = MERGE_TAXONOMY_BOLTZ2.out.csv_files
                .flatten()
                .map { f -> tuple(f.name.tokenize('.')[0], f) }

            all_csv = cached_csv.mix(computed_csv)
        }

        // --- Branch by predictor ---
        if (predictor == 'boltz2') {

            // Stage 4: Join MSA CSVs with bait-candidate pairs
            pairs_with_bait_msa = pairs_ch
                .map { pair_id, bait_id, bait_seq, candidate_id, candidate_seq ->
                    tuple(bait_id, pair_id, bait_seq, candidate_id, candidate_seq)
                }
                .combine(all_csv, by: 0)
                .map { bait_id, pair_id, bait_seq, candidate_id, candidate_seq, bait_csv ->
                    tuple(candidate_id, pair_id, bait_id, bait_seq, candidate_seq, bait_csv)
                }

            pairs_with_both_msa = pairs_with_bait_msa
                .combine(all_csv, by: 0)
                .map { candidate_id, pair_id, bait_id, bait_seq, candidate_seq, bait_csv, candidate_csv ->
                    tuple(pair_id, bait_id, bait_seq, candidate_id, candidate_seq, bait_csv, candidate_csv)
                }

            predict_ch = MAKE_BOLTZ_YAML(
                pairs_with_both_msa,
                ptms_path
            ).yaml_with_msa
                .map { pair_id, yaml, bait_csv, cand_csv -> tuple(pair_id, yaml, [bait_csv, cand_csv]) }

        } else {
            // predictor == 'alphafold3'

            // Stage 4: Join unpaired + paired A3M files with bait-candidate pairs
            pairs_with_bait_msa = pairs_ch
                .map { pair_id, bait_id, bait_seq, candidate_id, candidate_seq ->
                    tuple(bait_id, pair_id, bait_seq, candidate_id, candidate_seq)
                }
                .combine(all_a3m, by: 0)
                .combine(all_paired_msa, by: 0)
                .map { bait_id, pair_id, bait_seq, candidate_id, candidate_seq, bait_a3m, bait_paired_msa ->
                    tuple(candidate_id, pair_id, bait_id, bait_seq, candidate_seq, bait_a3m, bait_paired_msa)
                }

            pairs_with_both_msa = pairs_with_bait_msa
                .combine(all_a3m, by: 0)
                .combine(all_paired_msa, by: 0)
                .map { candidate_id, pair_id, bait_id, bait_seq, candidate_seq, bait_a3m, bait_paired_msa, candidate_a3m, candidate_paired_msa ->
                    tuple(pair_id, bait_id, bait_seq, candidate_id, candidate_seq, bait_a3m, candidate_a3m, bait_paired_msa, candidate_paired_msa)
                }

            predict_ch = MAKE_AF3_JSON(
                pairs_with_both_msa,
                ptms_path
            ).json_with_msa
                .map { pair_id, json_file, a3m_files -> tuple(pair_id, json_file, a3m_files) }
        }
    }


    // ---------------------------------------------------------------
    // Stage 5: Structure prediction (GPU-bound)
    // ---------------------------------------------------------------
    if (predictor == 'boltz2') {
        BOLTZ_PREDICT(predict_ch)
        pred_results = BOLTZ_PREDICT.out.results
    } else {
        AF3_PREDICT(predict_ch)
        pred_results = AF3_PREDICT.out.results
    }


    // ---------------------------------------------------------------
    // Stage 6: Parse metrics per pair (pDockQ + confidence)
    // ---------------------------------------------------------------
    if (predictor == 'boltz2') {
        PARSE_METRICS(pred_results)
        all_metrics = PARSE_METRICS.out.metrics.collect()
    } else {
        PARSE_METRICS_AF3(pred_results)
        all_metrics = PARSE_METRICS_AF3.out.metrics.collect()
    }


    // ---------------------------------------------------------------
    // Stage 7: Aggregate and rank all predictions
    // ---------------------------------------------------------------
    RANK_PREDICTIONS(all_metrics)
}


workflow.onComplete {

    // -------------------------------------------------------------------
    // Reconcile pairs enumerated against pairs actually ranked.
    //
    // The deliberate process-wide errorStrategy in nextflow.config is
    //   { task.attempt <= 2 ? 'retry' : 'ignore' }
    // on a several-hundred-candidate screen, one pair prediction that has
    // and OOM must not take the rest of the run down with it. The cost of
    // 'ignore' is that workflow.success stays true no matter how many
    // predictions were dropped.
    //
    // So count rows rather than trust the exit status. Both files carry a
    // header row, hence the -1. A null return means the file is absent,
    // which for the rankings file IS the diagnosis.
    // -------------------------------------------------------------------
    def dataRows = { pathStr ->
        def f = file(pathStr)
        f.exists() ? Math.max(0, f.readLines().count { it.trim() } - 1) : null
    }

    def n_pairs  = dataRows("${params.outdir}/pairs/pairs.tsv")
    def n_ranked = dataRows("${params.outdir}/rankings_${params.predictor}.tsv")

    // workflow.stats.ignoredCount counts TASKS swallowed by 'ignore', not
    // pairs: a pair that failed, retried twice and was then dropped
    // contributes one ignored task, and a pair can also be lost in
    // PARSE_METRICS rather than in prediction. Report it as a pointer to
    // the trace file, not as a pair count: n_pairs/n_ranked above are the
    // pair-level numbers.
    def n_ignored = workflow.stats.ignoredCount

    def status
    if (!workflow.success)
        // A real pipeline failure (not an ignored task): Nextflow has
        // already printed the cause above.
        status = 'FAILED'
    else if (n_pairs == null)
        // No pairs.tsv published, so there is nothing to reconcile against,
        // e.g. GENERATE_PAIRS itself was the thing that died.
        status = 'SUCCESS (no pair list found to reconcile against)'
    else if (n_ranked == null)
        status = "FAILED - no rankings written (0 of ${n_pairs} pairs ranked)"
    else if (n_ranked < n_pairs)
        status = "INCOMPLETE - ${n_ranked} of ${n_pairs} pairs ranked, ${n_pairs - n_ranked} dropped"
    else
        status = "SUCCESS - ${n_ranked} of ${n_pairs} pairs ranked"

    // Assembled line by line rather than as one GString so the optional
    // "Dropped" line can be omitted without leaving a blank row or
    // disturbing stripIndent()'s common-prefix calculation.
    def lines = []
    lines << '========================================================================='
    lines << 'OmegaPulldown complete'
    lines << '========================================================================='
    lines << "Predictor  : ${params.predictor}"
    lines << "Status     : ${status}"
    lines << "Duration   : ${workflow.duration}"
    if (n_ignored > 0)
        lines << "Dropped    : ${n_ignored} task(s) failed after retries and were ignored"
    lines << "Output     : ${params.outdir}"
    lines << "Structures : ${params.outdir}/structures_${params.predictor}/"
    lines << "Rankings   : ${params.outdir}/rankings_${params.predictor}.tsv"
    lines << '========================================================================='

    log.info '\n' + lines.collect { '    ' + it }.join('\n') + '\n'
}
