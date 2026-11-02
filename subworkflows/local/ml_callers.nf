// Machine-learning callers.
//  LYCEUM    pretrained (released weights); runs in every mode.
//  CNV-MEANN and dudeML must be trained on labelled data:
//    --mode simulate : trained here with leave-one-replicate-out, as in the paper (each CNV-positive sample is
//                      predicted by models that never saw its replicate; negative controls by all models), and the
//                      models are saved to <outdir>/models/
//    other modes     : applied with models given by --ml_models (e.g. <outdir>/models/ of a simulate run),
//                      as for Simulation 2 in the paper; skipped with a warning when no models are given.
// emits: results = [ meta, tool, config, 'calls', file ]

include { LYCEUM; CNVMEANN_FEATURES; CNVMEANN_TRAIN; CNVMEANN_PREDICT;
          DUDEML_FEATURES; DUDEML_TRAIN; DUDEML_PREDICT } from '../../modules/local/ml'

def kbl(b) { int x = b as int; "${(x / 1000) as int}kb".toString() }

workflow ML_CALLERS {
    take:
    ch_bam     // [ meta, bam, bai ]   (simulate: meta.rep = replicate, 0 for negative controls)
    ch_truth   // [ rep, truth ]       (simulate only; Channel.empty() otherwise)
    ch_ref     // [ fasta, fai ]
    do_train   // true in simulate mode

    main:
    callers = params.callers.toString().toLowerCase().tokenize(',')*.trim()
    NO_FILE = file("${projectDir}/assets/NO_FILE")
    reps    = (1..(params.sim_replicates as int)).toList()
    ch_out = Channel.empty()

    if ('lyceum' in callers) {
        chain = params.genome == 'GRCh37' ?
            file(params.liftover_chain ?: 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz') : NO_FILE
        LYCEUM(ch_bam, chain)
        ch_out = ch_out.mix(LYCEUM.out.calls.map { m, c, f -> [m, 'LYCEUM', c, 'calls', f] })
    }

    trained = callers.findAll { it in ['cnvmeann', 'dudeml'] }
    if (trained && !do_train && !params.ml_models) {
        log.warn "CNV-MEANN and dudeML need trained models: run --mode simulate first and pass its <outdir>/models as --ml_models. Skipping ${trained.join(' and ')}."
    }
    use_trained = do_train || params.ml_models

    // labelled input: CNV-positive simulated samples carry their truth set, everything else NO_FILE
    if (do_train) {
        ch_pos = ch_bam.filter { it[0].rep }.map { m, b, i -> [m.rep, m, b, i] }
                       .combine(ch_truth, by: 0).map { r, m, b, i, t -> [m, b, i, t] }
        ch_neg = ch_bam.filter { !it[0].rep }.map { m, b, i -> [m, b, i, NO_FILE] }
        ch_in  = ch_pos.mix(ch_neg)
    } else {
        ch_in  = ch_bam.map { m, b, i -> [m, b, i, NO_FILE] }
    }

    if ('cnvmeann' in callers && use_trained) {
        CNVMEANN_FEATURES(ch_in, ch_ref)
        F = CNVMEANN_FEATURES.out.features
        if (do_train) {
            // one network per held-out replicate (trained on the other replicates, all depths) + one on everything
            ch_tr = F.filter { it[0].rep }.map { m, d -> [m.rep, d] }.collect(flat: false).flatMap { L ->
                reps.collect { f -> ["fold${f}".toString(), L.findAll { it[0] != f }.collect { it[1] }] } +
                [["all", L.collect { it[1] }]] }
            CNVMEANN_TRAIN(ch_tr)
            ch_pred = F.map { m, d -> [(m.rep ? "fold${m.rep}" : "all").toString(), m, d] }
                       .combine(CNVMEANN_TRAIN.out.model, by: 0).map { k, m, d, f -> [m, d, f] }
        } else {
            mf = file("${params.ml_models}/CNV-MEANN/CNV-MEANN_model_all.joblib", checkIfExists: true)
            ch_pred = F.map { m, d -> [m, d, mf] }
        }
        CNVMEANN_PREDICT(ch_pred, kbl(params.cnvmeann_bin))
        ch_out = ch_out.mix(CNVMEANN_PREDICT.out.calls.map { m, c, f -> [m, 'CNV-MEANN', c, 'calls', f] })
    }

    if ('dudeml' in callers && use_trained) {
        DUDEML_FEATURES(ch_in)
        S = DUDEML_FEATURES.out.sample
        if (do_train) {
            // per depth: one random forest per held-out replicate, trained on the other replicates at that depth
            ch_tr = DUDEML_FEATURES.out.train.map { m, f -> [m.level, m.rep, f] }.groupTuple(by: 0)
                .flatMap { lvl, rs, fs ->
                    def pairs = [rs, fs].transpose()
                    reps.findAll { f -> pairs.any { it[0] != f } }
                        .collect { f -> [lvl, lvl.toString().replace('x', ''), f, pairs.findAll { it[0] != f }.collect { it[1] }] } }
            DUDEML_TRAIN(ch_tr)
            byLevel = DUDEML_TRAIN.out.model.groupTuple(by: 0)          // [ level, [folds], [[sav, sav2], ...] ]
            ch_pred = S.map { m, b, d -> [m.level, m, b, d] }.combine(byLevel, by: 0).map { l, m, b, d, folds, files ->
                def sel = m.rep ? files[folds.indexOf(m.rep)] : files
                [m, b, d, [sel].flatten()] }
        } else {
            md = file("${params.ml_models}/dudeML", checkIfExists: true)
            ch_pred = S.map { m, b, d -> [m, b, d, md] }
        }
        DUDEML_PREDICT(ch_pred, kbl(params.dudeml_bin))
        ch_out = ch_out.mix(DUDEML_PREDICT.out.calls.map { m, c, f -> [m, 'dudeML', c, 'calls', f] })
    }

    emit:
    results = ch_out
}
