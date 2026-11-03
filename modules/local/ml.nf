// Machine-learning callers: CNV-MEANN and dudeML (trained by the workflow) and LYCEUM (pretrained).

process LYCEUM {
    tag "${meta.id}"
    label 'cfcnv_lyceum'
    cpus 8
    memory { 32.GB * task.attempt }
    time { 12.h * task.attempt }

    input:
    tuple val(meta), path(bam), path(bai)
    path chain

    output:
    tuple val(meta), val('exon'), path("${meta.id}.LYCEUM.exon.segments.tsv"), emit: calls
    path "${meta.id}.LYCEUM.exon.segments.exonlevel.csv", emit: exons

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    run_lyceum.sh ${bam} ${params.genome} ${chain} ${task.cpus} ${meta.id}.LYCEUM.exon.segments.tsv
    """
    stub:
    """
    printf 'chrom\\tstart\\tend\\ttype\\tn_exons\\n1\\t10000000\\t10500000\\tloss\\t3\\n' > ${meta.id}.LYCEUM.exon.segments.tsv
    touch ${meta.id}.LYCEUM.exon.segments.exonlevel.csv
    """
}

process CNVMEANN_FEATURES {
    tag "${meta.id}"
    label 'cfcnv_ml'
    cpus 2
    memory { 24.GB * task.attempt }
    time { 12.h * task.attempt }

    input:
    tuple val(meta), path(bam), path(bai), path(truth)
    tuple path(fasta), path(fai)

    output:
    tuple val(meta), path("${meta.id}_meann"), emit: features

    script:
    def B = params.cnvmeann_bin as int
    def t = truth.name != 'NO_FILE' ? "--truth ${truth}" : ''
    """
    cnvmeann_features.py --bam ${bam} --ref ${fasta} --sample ${meta.id} --outdir raw --bin ${B} ${t}
    cnvmeann_normalize.py --raw raw --out ${meta.id}_meann --sample ${meta.id} --bin ${B}
    rm -rf raw
    """
    stub:
    """
    mkdir -p ${meta.id}_meann && touch ${meta.id}_meann/.done
    """
}

process CNVMEANN_TRAIN {
    tag "${key}"
    label 'cfcnv_ml'
    cpus 4
    memory { 48.GB * task.attempt }
    time { 8.h * task.attempt }

    input:
    tuple val(key), path(dirs, stageAs: 'train/*')

    output:
    tuple val(key), path("CNV-MEANN_model_${key}.joblib"), emit: model

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    cnvmeann_train.py --bin ${params.cnvmeann_bin} --out CNV-MEANN_model_${key}.joblib train/*
    """
    stub:
    """
    touch CNV-MEANN_model_${key}.joblib
    """
}

process CNVMEANN_PREDICT {
    tag "${meta.id}"
    label 'cfcnv_ml'
    label 'small'
    memory { 8.GB * task.attempt }

    input:
    tuple val(meta), path(features), path(model)
    val cfg

    output:
    tuple val(meta), val(cfg), path("${meta.id}.CNV-MEANN.${cfg}.events.tsv"), emit: calls

    script:
    """
    cnvmeann_predict.py --model ${model} --features ${features} --out ${meta.id}.CNV-MEANN.${cfg}.events.tsv
    """
    stub:
    """
    printf 'chrom\\tstart\\tend\\ttype\\n1\\t10000000\\t11000000\\tloss\\n' > ${meta.id}.CNV-MEANN.${cfg}.events.tsv
    """
}

process DUDEML_FEATURES {
    tag "${meta.id}"
    label 'cfcnv_ml'
    cpus 2
    memory { 16.GB * task.attempt }
    time { 12.h * task.attempt }

    input:
    tuple val(meta), path(bam), path(bai), path(truth)

    output:
    tuple val(meta), path("${meta.id}.sample.bed"), path("${meta.id}.depth.txt"), emit: sample
    tuple val(meta), path("${meta.id}.train.bed"), emit: train, optional: true

    script:
    def t = truth.name != 'NO_FILE' ? truth : 'NONE'
    """
    dudeml_features.sh ${bam} ${meta.id} ${params.dudeml_bin} ${t}
    estimate_depth.sh ${bam} > ${meta.id}.depth.txt
    """
    stub:
    def tr = truth.name != 'NO_FILE' ? "touch ${meta.id}.train.bed" : ''
    """
    touch ${meta.id}.sample.bed; echo 1 > ${meta.id}.depth.txt; ${tr}
    """
}

process DUDEML_TRAIN {
    tag "${level} fold ${fold}"
    label 'cfcnv_ml'
    cpus 4
    memory { 32.GB * task.attempt }
    time { 8.h * task.attempt }

    input:
    tuple val(level), val(d), val(fold), path(beds, stageAs: 'train/*')

    output:
    tuple val(level), val(fold), path("depth_${d}/fold${fold}.sav*"), emit: model

    script:
    """
    mkdir -p depth_${d}
    dudeml_train.sh depth_${d}/fold${fold}.sav ${2026 + (fold as int)} train/*
    """
    stub:
    """
    mkdir -p depth_${d} && touch depth_${d}/fold${fold}.sav depth_${d}/fold${fold}.sav2
    """
}

process DUDEML_PREDICT {
    tag "${meta.id}"
    label 'cfcnv_ml'
    label 'small'
    memory { 8.GB * task.attempt }

    input:
    tuple val(meta), path(sample_bed), path(depth), path(models, stageAs: 'models/*')
    val cfg

    output:
    tuple val(meta), val(cfg), path("${meta.id}.dudeML.${cfg}.predict.bed.gz"), emit: calls

    script:
    """
    # models/ holds model pairs (fold*.sav + fold*.sav2), or a model folder with depth_<x>/ sub-folders:
    # then the folder trained at the depth closest to this sample's depth is used (all its fold models vote)
    M=models; [ -d models/dudeML ] && M=models/dudeML
    if ls -d \$M/depth_* >/dev/null 2>&1; then
      D=\$(cat ${depth})
      M=\$(for x in \$M/depth_*; do echo "\$x \${x##*depth_}"; done | awk -v d="\$D" '{e=\$2-d; if(e<0)e=-e; print e, \$1}' | sort -g | head -n1 | cut -d' ' -f2)
      echo "sample depth \${D}x -> models in \$M"
    fi
    D=\${DUDEML_DIR:-\$(get_dudeml.sh)}
    python \$D/dudeML.py predict -i ${sample_bed} -t \$(readlink -f \$M)/ -o ${meta.id}.dudeML.${cfg}.predict.bed
    gzip -f ${meta.id}.dudeML.${cfg}.predict.bed
    """
    stub:
    """
    printf '1\\t10000001\\t11000000\\tS\\t0.5\\tDel\\t0.9\\t1\\t0.9\\n' | gzip > ${meta.id}.dudeML.${cfg}.predict.bed.gz
    """
}
