// Run the selected callers and tumour-fraction estimators on a channel of BAMs.
// emits: results = [ meta, tool, config, kind (calls|tf), file ]

include { QDNASEQ; ACE; ICHORCNA; CNVKIT_REFERENCE; CNVKIT; MOSDEPTH_BINS; SPECTRE;
          BICSEQ2_RESOURCES; BICSEQ2_UNIQUE; BICSEQ2_NORMSEG; FRAGLE } from '../../modules/local/callers'

def kbLabel(b) { int x = b as int; x >= 1000 ? "${(x / 1000) as int}kb" : "${x}bp" }

workflow CALL_CNAS {
    take:
    ch_bam      // [ meta, bam, bai ]
    ch_ref      // [ fasta, fai ]  (value channel)
    fai         // path to .fai (value)

    main:
    def callers = params.callers.toString().toLowerCase().tokenize(',')*.trim()
    def tfm     = params.tf_methods.toString().toLowerCase().tokenize(',')*.trim()
    ch_out = Channel.empty()

    if ('qdnaseq' in callers || 'ace' in tfm) {
        def bins = params.qdnaseq_bins_kb.toString().tokenize(',')*.trim()
        QDNASEQ(ch_bam.combine(Channel.fromList(bins)))
        if ('qdnaseq' in callers) ch_out = ch_out.mix(QDNASEQ.out.calls.map { m, c, f -> [m, 'QDNAseq', c, 'calls', f] })
        if ('ace' in tfm) {
            ACE(QDNASEQ.out.rds)
            ch_out = ch_out.mix(ACE.out.tf.map { m, c, f -> [m, 'ACE', c, 'tf', f] })
        }
    }
    if ('ichorcna' in callers || 'ichorcna' in tfm) {
        def modes = params.ichor_modes.toString().tokenize(',')*.trim()
        def W  = params.ichor_window_kb as int
        def wl = W >= 1000 ? "${(W / 1000) as int}Mb" : "${W}kb"
        def ml = modes.collect { m -> [m, "${wl} ${m == 'lowTF' ? 'low-TF' : 'defaults'}".toString(), "${wl}_${m}".toString()] }
        ICHORCNA(ch_bam.combine(Channel.fromList(ml)), fai)
        if ('ichorcna' in callers) ch_out = ch_out.mix(ICHORCNA.out.calls.map { m, c, f -> [m, 'ichorCNA', c, 'calls', f] })
        if ('ichorcna' in tfm)     ch_out = ch_out.mix(ICHORCNA.out.tf.map    { m, c, f -> [m, 'ichorCNA', c, 'tf', f] })
    }
    if ('cnvkit' in callers) {
        if (params.cnvkit_normals) {
            ch_normals = Channel.fromPath(params.cnvkit_normals.toString().tokenize(',')*.trim(), checkIfExists: true)
                .flatMap { b -> [b, file("${b}.bai", checkIfExists: true)] }.collect()
            CNVKIT_REFERENCE(ch_normals, ch_ref); ch_cnn = CNVKIT_REFERENCE.out.cnn
        } else {
            ch_cnn = Channel.value(file("${projectDir}/assets/NO_FILE"))
        }
        CNVKIT(ch_bam, ch_cnn, ch_ref, kbLabel(params.cnvkit_bin))
        ch_out = ch_out.mix(CNVKIT.out.calls.map { m, c, f -> [m, 'CNVkit', c, 'calls', f] })
    }
    if ('spectre' in callers) {
        bl = file(params.blacklist ?: "${projectDir}/assets/NO_FILE", checkIfExists: true)
        MOSDEPTH_BINS(ch_bam)
        SPECTRE(MOSDEPTH_BINS.out.bed, ch_ref, bl, kbLabel(params.spectre_bin))
        ch_out = ch_out.mix(SPECTRE.out.calls.map { m, c, f -> [m, 'Spectre', c, 'calls', f] })
    }
    if ('bicseq2' in callers && !params.bicseq2_container && !workflow.stubRun) {
        log.warn "BIC-seq2 skipped: cfCNV does not distribute a BIC-seq2 image (academic-use licence). Build one as described in assets/hpc/README_bicseq2.md and pass --bicseq2_container."
    }
    if ('bicseq2' in callers && (params.bicseq2_container || workflow.stubRun)) {
        if (params.bicseq2_resources) {
            ch_bres = Channel.value(file(params.bicseq2_resources, checkIfExists: true))
        } else {
            log.warn "No --bicseq2_resources given: building BIC-seq2 mappability resources first (slow, ~1-2 days, done once)."
            BICSEQ2_RESOURCES(ch_ref); ch_bres = BICSEQ2_RESOURCES.out.dir
        }
        BICSEQ2_UNIQUE(ch_bam, fai)
        BICSEQ2_NORMSEG(BICSEQ2_UNIQUE.out.unique, ch_bres, fai, kbLabel(params.bicseq2_bin))
        ch_out = ch_out.mix(BICSEQ2_NORMSEG.out.calls.map { m, c, f -> [m, 'BIC-seq2', c, 'calls', f] })
    }
    if ('fragle' in tfm) {
        FRAGLE(ch_bam)
        ch_out = ch_out.mix(FRAGLE.out.tf.map { m, c, f -> [m, 'Fragle', c, 'tf', f] })
    }

    emit:
    results = ch_out
}
