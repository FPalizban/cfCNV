// Harmonisation, consensus, annotation and scoring

process COLLECT_CALLS {
    label 'cfcnv_tools'
    label 'small'
    memory { 8.GB * task.attempt }

    input:
    path manifest
    path files, stageAs: 'inputs/*'
    path fai

    output:
    path 'all_calls.tsv',        emit: calls
    path 'call_burden.tsv',      emit: burden
    path 'tumour_fraction.tsv',  emit: tf

    script:
    """
    # manifest paths are file names; point them at the staged copies
    awk 'BEGIN{FS=OFS="\\t"} NR==1{print; next} {\$9="inputs/"\$9; print}' ${manifest} > manifest.staged.tsv
    collect_calls.py --manifest manifest.staged.tsv --fai ${fai} --bic-threshold ${params.bicseq2_log2_threshold} --outdir .
    """
    stub:
    """
    awk 'BEGIN{FS=OFS="\\t"} NR==1{print; next} {\$9="inputs/"\$9; print}' ${manifest} > manifest.staged.tsv
    collect_calls.py --manifest manifest.staged.tsv --fai ${fai} --outdir .
    """
}

process CONSENSUS {
    label 'cfcnv_tools'
    label 'small'
    memory { 8.GB * task.attempt }

    input:
    path calls
    path burden
    path fai
    path samplesheet
    path tcga
    val primary

    output:
    path 'caller_agreement.tsv',  emit: agreement
    path 'consensus_summary.tsv', emit: summary
    path 'consensus_calls.tsv',   emit: consensus
    path 'depth_robustness.tsv',  emit: robustness
    path 'window_states.tsv.gz',  emit: states, optional: true

    script:
    def t = tcga.name != 'NO_FILE' ? "--tcga ${tcga}" : ''
    """
    consensus.py --calls ${calls} --fai ${fai} --primary "${primary}" --samples ${samplesheet} ${t} \\
        --window ${params.consensus_window} --min-callers ${params.consensus_min} --outdir .
    """
    stub:
    def t = tcga.name != 'NO_FILE' ? "--tcga ${tcga}" : ''
    """
    consensus.py --calls ${calls} --fai ${fai} --primary "${primary}" --samples ${samplesheet} ${t} \\
        --window ${params.consensus_window} --min-callers ${params.consensus_min} --outdir .
    """
}

process ANNOTATE {
    label 'cfcnv_tools'
    label 'medium'

    input:
    path calls
    path refs
    path chain
    val primary

    output:
    path 'annotation/*', emit: tables

    script:
    def lift = chain.name != 'NO_FILE' ? "--lift ${chain}" : ''
    """
    annotate_calls.py ${calls} ${refs} annotation ${lift} --configs "${primary}"
    """
    stub:
    """
    mkdir -p annotation && touch annotation/calls_annotated.tsv
    """
}

process SCORE {
    label 'cfcnv_tools'
    label 'small'
    memory { 8.GB * task.attempt }

    input:
    path calls
    path burden
    path truth_map
    path truth_files, stageAs: 'truth/*'
    path fai

    output:
    path 'event_scores.tsv',    emit: scores
    path 'truth_detection.tsv', emit: detection

    script:
    """
    awk 'BEGIN{FS=OFS="\\t"} {if(\$2!="" && \$2!="NA") \$2="truth/"\$2; print}' ${truth_map} > truth_map.staged.tsv
    score_calls.py --calls ${calls} --truth-map truth_map.staged.tsv --fai ${fai} --overlap ${params.score_overlap} --sample user --outdir .
    """
    stub:
    """
    awk 'BEGIN{FS=OFS="\\t"} {if(\$2!="" && \$2!="NA") \$2="truth/"\$2; print}' ${truth_map} > truth_map.staged.tsv
    score_calls.py --calls ${calls} --truth-map truth_map.staged.tsv --fai ${fai} --overlap ${params.score_overlap} --sample user --outdir .
    """
}
