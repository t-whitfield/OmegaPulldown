/*
 * =========================================================================
 *  OmegaPulldown: Protein-protein interaction screen using Boltz-2
 *                 or AlphaFold 3
 * =========================================================================
 *
 *  Analogous to AlphaPulldown (see D. Yu et al. Bioinformatics 39, 2023,
 *  btac749) but using Boltz-2 (cite S. Passaro et al. bioRxiv, 2025,
 *  DOI:10.1101/2025.06.14.659707) or AlphaFold 3 (cite J. Abramson et al.
 *  Nature 630, 2024, 493-500) for structure prediction.
 *
 *  Input:  two multi-sequence FASTA files (baits and candidates)
 *  Output: ranked table of predicted complex structures with confidence
 *          metrics (ipTM, pDockQ, pLDDT, pTM), plus per-pair structure
 *          subdirectories under results/structures/
 *
 *  Pipeline stages:
 *    1. GENERATE_PAIRS    - enumerate bait x candidate combinations
 *    2. MMSEQS_SEARCH     - compute MSAs per unique protein (CPU)
 *                           [skipped if pre-computed A3Ms exist]
 *    3a. A3M_TO_BOLTZ_CSV - convert a3m to Boltz-2 paired CSV format
 *    3b. (AF3 path)       - A3M files used directly via unpairedMsaPath
 *    4. MAKE_BOLTZ_YAML / MAKE_AF3_JSON - assemble input per pair
 *    5. BOLTZ_PREDICT / AF3_PREDICT     - structure prediction (GPU)
 *    6. PARSE_METRICS     - extract confidence + pDockQ per pair
 *    7. RANK_PREDICTIONS  - aggregate and rank all pairs
 *
 *  Predictor selection:
 *    --predictor boltz2       (default)
 *    --predictor alphafold3
 *
 *  MSA caching:
 *    If A3M files exist in params.alignments_dir/a3m/<prot_id>.a3m,
 *    MMSEQS_SEARCH is skipped for those proteins.
 *
 *  When params.use_msa_server = true, stages 2-3 are skipped and
 *  Boltz-2 fetches MSAs from the MMseqs2 server at prediction time.
 *  (Only valid with predictor = boltz2.)
 *
 *  Post-translational modifications (optional):
 *    Place bait_ptms.csv and/or candidate_ptms.csv in the same directory
 *    as the FASTA files, or specify paths explicitly via --bait_ptms and
 *    --candidate_ptms.  CSV format:
 *      protein_id,position,ccd
 *      BAIT1_GFP,18,SEP
 *
 *  Usage:
 *    # Boltz-2 (default) on a Slurm cluster with GPUs
 *    nextflow run main.nf -profile slurm \
 *      --baits_fasta inputs/sequences/baits.fasta \
 *      --candidates_fasta inputs/sequences/candidates.fasta \
 *      --mmseqs_db /path/to/uniref90
 *
 *    # AlphaFold 3
 *    nextflow run main.nf -profile slurm \
 *      --predictor alphafold3 \
 *      --baits_fasta inputs/sequences/baits.fasta \
 *      --candidates_fasta inputs/sequences/candidates.fasta
 *
 *    # Skip MSA computation (use pre-computed alignments)
 *    nextflow run main.nf -profile slurm \
 *      --predictor alphafold3 \
 *      --alignments_dir results/alignments
 * =========================================================================
 */

nextflow.enable.dsl = 2

// --- Module imports ---
include { MMSEQS_SEARCH; A3M_TO_BOLTZ_CSV }  from './modules/search_msa'
include { GENERATE_PAIRS; MAKE_BOLTZ_YAML; MAKE_BOLTZ_YAML_NO_MSA }  from './modules/prepare_input'
include { BOLTZ_PREDICT }  from './modules/predict'
include { PARSE_METRICS; RANK_PREDICTIONS }  from './modules/rank'
include { MAKE_AF3_JSON; AF3_PREDICT }  from './modules/predict_af3'
include { PARSE_METRICS_AF3 }  from './modules/rank_af3'


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
 * Resolve an optional PTM CSV path.
 *
 * Priority:
 *   1. Explicit parameter (--bait_ptms / --candidate_ptms)
 *   2. Convention: <conventionName> in the same directory as the FASTA
 *   3. No PTMs (returns 'no_file' sentinel)
 *
 * Returns a path string suitable for passing to make_boltz_yaml.py.
 */
def resolvePtmsCsv(explicitParam, fastaPath, conventionName) {
    // 1. Explicit parameter
    if (explicitParam != null) {
        def f = file(explicitParam)
        if (f.exists()) {
            log.info "PTM CSV (explicit): ${f}"
            return f
        } else {
            log.warn "PTM CSV specified but not found: ${explicitParam}"
            return 'no_file'
        }
    }
    // 2. Convention: role-specific CSV alongside the FASTA
    if (fastaPath == null) {
        return 'no_file'
    }
    def fastaDir = file(fastaPath).parent
    def conventionPath = fastaDir.resolve(conventionName)
    if (conventionPath.exists()) {
        log.info "PTM CSV (auto-detected): ${conventionPath}"
        return conventionPath
    }
    // 3. No PTMs
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
        error "use_msa_server is not supported with predictor=alphafold3. " +
              "Provide pre-computed MSAs or use local mmseqs2 search."
    }

    log.info """
    =========================================================================
    boltzPulldown
    =========================================================================
    Predictor    : ${predictor}
    Baits        : ${params.baits_fasta}
    Candidates   : ${params.candidates_fasta}
    Alignments   : ${params.alignments_dir}
    MSA server   : ${params.use_msa_server}
    Output       : ${params.outdir}
    =========================================================================
    """.stripIndent()

    // ---------------------------------------------------------------
    // Resolve optional PTM CSV files
    // ---------------------------------------------------------------
    bait_ptms_path      = resolvePtmsCsv(params.bait_ptms, params.baits_fasta, 'bait_ptms.csv')
    candidate_ptms_path = resolvePtmsCsv(params.candidate_ptms, params.candidates_fasta, 'candidate_ptms.csv')

    // ---------------------------------------------------------------
    // Stage 1: Generate all bait x candidate pairs
    // ---------------------------------------------------------------
    baits_ch      = Channel.fromPath(params.baits_fasta)
    candidates_ch = Channel.fromPath(params.candidates_fasta)

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
            bait_ptms_path,
            candidate_ptms_path
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

        // --- MSA caching: split into cached vs. need-to-compute ---
        unique_proteins_ch.branch {
            cached:  file("${params.alignments_dir}/a3m/${it[0]}.a3m").exists()
            compute: true
        }.set { split_ch }

        // Use cached A3M files directly
        cached_a3m = split_ch.cached.map { prot_id, seq ->
            tuple(prot_id, file("${params.alignments_dir}/a3m/${prot_id}.a3m"))
        }

        // Stage 2: MSA search for proteins without cached A3Ms (CPU-bound)
        MMSEQS_SEARCH(split_ch.compute)

        // Merge cached + newly computed A3M streams
        all_a3m = cached_a3m.mix(MMSEQS_SEARCH.out.a3m)

        // --- Branch by predictor ---
        if (predictor == 'boltz2') {

            // Stage 3: Convert a3m -> Boltz-2 CSV with pairing keys
            // Also check for cached CSV files
            all_a3m.branch {
                cached:  file("${params.alignments_dir}/csv/${it[0]}.csv").exists()
                compute: true
            }.set { csv_split_ch }

            cached_csv = csv_split_ch.cached.map { prot_id, a3m ->
                tuple(prot_id, file("${params.alignments_dir}/csv/${prot_id}.csv"))
            }

            A3M_TO_BOLTZ_CSV(csv_split_ch.compute)

            all_csv = cached_csv.mix(A3M_TO_BOLTZ_CSV.out.csv)

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
                bait_ptms_path,
                candidate_ptms_path
            ).yaml_with_msa
                .map { pair_id, yaml, bait_csv, cand_csv -> tuple(pair_id, yaml, [bait_csv, cand_csv]) }

        } else {
            // predictor == 'alphafold3'

            // Stage 4: Join A3M files with bait-candidate pairs
            pairs_with_bait_msa = pairs_ch
                .map { pair_id, bait_id, bait_seq, candidate_id, candidate_seq ->
                    tuple(bait_id, pair_id, bait_seq, candidate_id, candidate_seq)
                }
                .combine(all_a3m, by: 0)
                .map { bait_id, pair_id, bait_seq, candidate_id, candidate_seq, bait_a3m ->
                    tuple(candidate_id, pair_id, bait_id, bait_seq, candidate_seq, bait_a3m)
                }

            pairs_with_both_msa = pairs_with_bait_msa
                .combine(all_a3m, by: 0)
                .map { candidate_id, pair_id, bait_id, bait_seq, candidate_seq, bait_a3m, candidate_a3m ->
                    tuple(pair_id, bait_id, bait_seq, candidate_id, candidate_seq, bait_a3m, candidate_a3m)
                }

            predict_ch = MAKE_AF3_JSON(
                pairs_with_both_msa,
                bait_ptms_path,
                candidate_ptms_path
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
    log.info """
    =========================================================================
    boltzPulldown complete
    =========================================================================
    Predictor  : ${params.predictor}
    Status     : ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration   : ${workflow.duration}
    Output     : ${params.outdir}
    Structures : ${params.outdir}/structures_${params.predictor}/
    Rankings   : ${params.outdir}/rankings_${params.predictor}.tsv
    =========================================================================
    """.stripIndent()
}
