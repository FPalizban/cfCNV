// BAM preparation: down-sampling, Simulation 1 spike-in, Simulation 2 mixtures, truth sets

process DOWNSAMPLE {
    tag "${meta.id} -> ${target}x"
    label 'cfcnv_tools'
    label 'small'
    cpus 2

    input:
    tuple val(meta), path(bam), path(bai), val(target)
    path fai

    output:
    tuple val(meta), val(target), path("${meta.sample}_${target}x.bam"), path("${meta.sample}_${target}x.bam.bai"), emit: bam
    path "${meta.sample}_${target}x.bam.depth.txt", emit: depth

    script:
    def seed = 11 + (meta.seed ?: 0)
    """
    downsample_bam.sh ${bam} ${target} ${params.genome_size_bp ?: fai} ${meta.sample}_${target}x.bam ${seed} ${task.cpus}
    """
    stub:
    """
    touch ${meta.sample}_${target}x.bam ${meta.sample}_${target}x.bam.bai
    printf 'source_depth\\ttarget_depth\\tfraction\\n10\\t${target}\\t0.1\\n' > ${meta.sample}_${target}x.bam.depth.txt
    """
}

process MAKE_TRUTH {
    tag "rep${rep}"
    label 'cfcnv_tools'
    label 'small'

    input:
    val rep
    path sizes

    output:
    tuple val(rep), path("rep${rep}.truth.tsv"), emit: truth

    script:
    """
    make_truth.py --chrom-sizes ${sizes} --replicate ${rep} --seed-base 1000 --output rep${rep}.truth.tsv
    """
    stub:
    """
    printf '#chrom\\tstart\\tend\\tcnv_type\\n1\\t10000000\\t11000000\\tloss\\n' > rep${rep}.truth.tsv
    """
}

process SPIKE_CNVS {
    tag "${meta.sample}_rep${rep}"
    label 'cfcnv_tools'
    label 'medium'

    input:
    tuple val(meta), path(bam), path(bai), val(rep), path(truth)
    path fai

    output:
    tuple val(meta), val(rep), path("${meta.sample}_rep${rep}.bam"), path("${meta.sample}_rep${rep}.bam.bai"), emit: bam

    script:
    """
    spike_cnvs.sh ${bam} ${truth} ${fai} ${meta.sample}_rep${rep}.bam ${1000 + (rep as int)} ${task.cpus}
    """
    stub:
    """
    touch ${meta.sample}_rep${rep}.bam ${meta.sample}_rep${rep}.bam.bai
    """
}

process MIX_TUMOUR_NORMAL {
    tag "TF${fraction} ${depth}x"
    label 'cfcnv_tools'
    label 'small'
    cpus 2

    input:
    tuple val(depth), val(fraction), val(seed)
    tuple path(tbam), path(tbai), path(nbam), path(nbai)

    output:
    tuple val(depth), val(fraction), path("mix_tf${fraction}_${depth}x.bam"), path("mix_tf${fraction}_${depth}x.bam.bai"), emit: bam

    script:
    """
    mix_tumour_normal.sh ${tbam} ${params.tumour_depth} ${nbam} ${params.normal_depth} ${depth} ${fraction} ${seed} \\
        mix_tf${fraction}_${depth}x.bam ${task.cpus}
    """
    stub:
    """
    touch mix_tf${fraction}_${depth}x.bam mix_tf${fraction}_${depth}x.bam.bai
    """
}
