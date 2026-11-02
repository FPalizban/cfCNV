#!/usr/bin/env nextflow
/*
 * cfCNV: benchmark of copy-number callers for low-pass cell-free DNA WGS
 * https://github.com/FPalizban/cfCNV
 *
 *   --mode call      run callers (+ optional down-sampling), harmonise, consensus, annotation   [default]
 *   --mode simulate  Simulation 1: spike events into healthy cfDNA BAMs, down-sample, call, train ML callers, score
 *   --mode mix       Simulation 2: tumour-normal mixtures at set depths and tumour fractions, call
 *   --mode score     score one call set (TSV/BED) against a truth set
 */
nextflow.enable.dsl = 2

include { DOWNSAMPLE; MAKE_TRUTH; SPIKE_CNVS; MIX_TUMOUR_NORMAL } from './modules/local/preprocess'
include { COLLECT_CALLS; CONSENSUS; ANNOTATE; SCORE }             from './modules/local/summarise'
include { CALL_CNAS }                                              from './subworkflows/local/call_cnas'
include { ML_CALLERS }                                             from './subworkflows/local/ml_callers'

def NO_FILE = file("${projectDir}/assets/NO_FILE")

// primary configuration per caller (used for consensus/annotation) = the first setting of each caller
def primaryConfigs() {
    def c = params.callers.toString().toLowerCase().tokenize(',')*.trim()
    def kb = { int b -> b >= 1000 ? "${(b / 1000) as int}kb" : "${b}bp" }
    def W = params.ichor_window_kb as int
    def wl = W >= 1000 ? "${(W / 1000) as int}Mb" : "${W}kb"
    def m = params.ichor_modes.toString().tokenize(',')*.trim()[0]
    def p = []
    if ('qdnaseq'  in c) p << "QDNAseq:${params.qdnaseq_bins_kb.toString().tokenize(',')[0].trim()}kb"
    if ('ichorcna' in c) p << "ichorCNA:${wl} ${m == 'lowTF' ? 'low-TF' : 'defaults'}"
    if ('bicseq2'  in c) p << "BIC-seq2:${kb(params.bicseq2_bin as int)}"
    if ('cnvkit'   in c) p << "CNVkit:${kb(params.cnvkit_bin as int)}"
    if ('spectre'  in c) p << "Spectre:${kb(params.spectre_bin as int)}"
    return p.join(',')
}

def readSamplesheet(path) {
    Channel.fromPath(path, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { r ->
            if (!r.sample || !r.bam) error "samplesheet needs the columns 'sample' and 'bam' (got: ${r.keySet()})"
            // relative paths: launch directory first, then the pipeline folder (used by -profile test)
            def res = { String p -> def f = file(p); (f.exists() || p.startsWith('/')) ? f : file("${projectDir}/${p}") }
            def bam = res(r.bam); if (!bam.exists()) error "BAM not found: ${r.bam}"
            def bai = r.bai ? res(r.bai) : (file("${bam}.bai").exists() ? file("${bam}.bai") : file(bam.toString().replaceAll(/\.bam$/, '.bai')))
            if (!bai.exists()) error "BAM index not found for ${r.bam} (run: samtools index ${r.bam})"
            def meta = [id: r.sample, sample: r.sample, level: 'full', group: r.group ?: '',
                        expected_tf: r.expected_tf ?: '', seed: 0]
            [meta, bam, bai]
        }
}

// write one manifest line per caller output and gather the files
workflow COLLECT {
    take:
    ch_results   // [ meta, tool, config, kind, file ]
    fai

    main:
    ch_manifest = ch_results
        .map { m, tool, cfg, kind, f -> [m.id, m.sample, m.level, m.group, m.expected_tf, tool, cfg, kind, f.name].join('\t') }
        .collectFile(name: 'manifest.tsv', newLine: true, sort: true,
                     seed: ['sample_id', 'base', 'level', 'group', 'expected_tf', 'tool', 'config', 'kind', 'file'].join('\t'))
    COLLECT_CALLS(ch_manifest, ch_results.map { it[4] }.collect(), fai)

    emit:
    calls  = COLLECT_CALLS.out.calls
    burden = COLLECT_CALLS.out.burden
    tf     = COLLECT_CALLS.out.tf
}

workflow {
    log.info """
    ===============================================
     cfCNV v${workflow.manifest.version}   mode: ${params.mode}
     callers : ${params.callers}
     TF      : ${params.tf_methods}
     genome  : ${params.genome}   fasta: ${params.fasta}
     outdir  : ${params.outdir}
    ===============================================
    """.stripIndent()

    if (params.mode == 'score') {
        if (!params.calls || !params.truth) error "--mode score needs --calls and --truth"
        if (!params.fasta) error "--fasta (or a .fai via --fasta) is needed for the genome size"
        fai = file("${params.fasta}.fai", checkIfExists: true)
        tm = Channel.of("user\t${file(params.truth).name}").collectFile(name: 'truth_map.tsv', newLine: true)
        SCORE(file(params.calls, checkIfExists: true), NO_FILE, tm, file(params.truth, checkIfExists: true), fai)
        return
    }

    if (!params.fasta) error "--fasta is required"
    fasta  = file(params.fasta, checkIfExists: true)
    fai    = file("${params.fasta}.fai", checkIfExists: true)
    ch_ref = Channel.value([fasta, fai])
    primary = primaryConfigs()
    targets = (params.downsample && !(params.downsample instanceof Boolean)) ? params.downsample.toString().tokenize(',')*.trim() : []

    if (params.mode == 'call') {
        if (!params.input) error "--input samplesheet.csv is required"
        ch_native = readSamplesheet(params.input)
        ch_bam = ch_native
        if (targets) {
            DOWNSAMPLE(ch_native.combine(Channel.fromList(targets)).map { m, b, i, t -> [m, b, i, t] }, fai)
            ch_bam = ch_native.mix(DOWNSAMPLE.out.bam.map { m, t, b, i ->
                [m + [id: "${m.sample}_${t}x", level: "${t}x"], b, i] })
        }
        CALL_CNAS(ch_bam, ch_ref, fai)
        ML_CALLERS(ch_bam, Channel.empty(), ch_ref, false)
        COLLECT(CALL_CNAS.out.results.mix(ML_CALLERS.out.results), fai); res = COLLECT.out
        if (primary.tokenize(',').size() >= 2) {
            tcga = params.annotation_refs ? file("${params.annotation_refs}/tcga_window_freq_1Mb.tsv") : NO_FILE
            if (params.genome != 'GRCh37' && tcga != NO_FILE) {
                log.warn "TCGA concordance uses GRCh37 windows; it is skipped for ${params.genome} calls (see annotation output instead)."
                tcga = NO_FILE
            }
            CONSENSUS(res.calls, res.burden, fai, file(params.input), tcga, primary)
        }
        if (params.annotation_refs) {
            chain = params.genome == 'GRCh37' ? NO_FILE :
                    file(params.liftover_chain ?: "${params.annotation_refs}/hg38ToHg19.over.chain.gz", checkIfExists: true)
            ANNOTATE(res.calls, file(params.annotation_refs, checkIfExists: true), chain, primary)
        }
    }

    else if (params.mode == 'simulate') {
        if (!params.input) error "--input samplesheet.csv with the healthy background BAMs is required"
        ch_bg = readSamplesheet(params.input)
        reps  = (1..(params.sim_replicates as int)).toList()
        depths = params.sim_depths.toString().tokenize(',')*.trim()
        if (params.truth) {
            ch_truth = Channel.fromList(reps).map { r -> [r, file("${params.truth}/rep${r}.truth.tsv", checkIfExists: true)] }
        } else {
            MAKE_TRUTH(Channel.fromList(reps), fai); ch_truth = MAKE_TRUTH.out.truth
        }
        SPIKE_CNVS(ch_bg.combine(ch_truth), fai)
        ch_pos = SPIKE_CNVS.out.bam.map { m, r, b, i -> [m + [sample: "${m.sample}_rep${r}", rep: r, seed: (r as int) * 100], b, i] }
        ch_neg = ch_bg.map { m, b, i -> [m + [rep: 0, seed: 5000], b, i] }
        DOWNSAMPLE(ch_pos.mix(ch_neg).combine(Channel.fromList(depths)), fai)
        ch_bam = DOWNSAMPLE.out.bam.map { m, t, b, i -> [m + [id: "${m.sample}_${t}x", level: "${t}x", group: m.rep ? 'positive' : 'negative'], b, i] }
        CALL_CNAS(ch_bam, ch_ref, fai)
        ML_CALLERS(ch_bam, ch_truth, ch_ref, true)
        COLLECT(CALL_CNAS.out.results.mix(ML_CALLERS.out.results), fai); res = COLLECT.out
        // truth map: sample_id -> truth file (empty for negative controls)
        ch_tm = ch_bam.map { m, b, i -> m.rep ? "${m.id}\trep${m.rep}.truth.tsv" : "${m.id}\t" }
                     .collectFile(name: 'truth_map.tsv', newLine: true, sort: true)
        SCORE(res.calls, res.burden, ch_tm, ch_truth.map { it[1] }.collect(), fai)
    }

    else if (params.mode == 'mix') {
        if (!(params.tumour_bam && params.normal_bam && params.tumour_depth && params.normal_depth))
            error "--mode mix needs --tumour_bam, --normal_bam, --tumour_depth and --normal_depth"
        def tb = file(params.tumour_bam, checkIfExists: true); def nb = file(params.normal_bam, checkIfExists: true)
        pools = Channel.value([tb, file("${tb}.bai", checkIfExists: true), nb, file("${nb}.bai", checkIfExists: true)])
        def fr = params.mix_fractions.toString().tokenize(',')*.trim()
        def dp = params.mix_depths.toString().tokenize(',')*.trim()
        def grid = []; int s = 1000
        dp.each { d -> fr.each { f -> grid << [d, f, s++] } }
        MIX_TUMOUR_NORMAL(Channel.fromList(grid), pools)
        ch_bam = MIX_TUMOUR_NORMAL.out.bam.map { d, f, b, i ->
            def id = "mix_tf${f}_${d}x"
            [[id: id, sample: id, level: "${d}x", group: 'mixture', expected_tf: f, seed: 0], b, i] }
        CALL_CNAS(ch_bam, ch_ref, fai)
        ML_CALLERS(ch_bam, Channel.empty(), ch_ref, false)
        COLLECT(CALL_CNAS.out.results.mix(ML_CALLERS.out.results), fai)
    }
    else error "unknown --mode ${params.mode} (call | simulate | mix | score)"
}

workflow.onComplete {
    log.info(workflow.success ? "\ncfCNV finished. Results in ${params.outdir}/" : "\ncfCNV failed - see .nextflow.log and the work directory of the failed task")
}
