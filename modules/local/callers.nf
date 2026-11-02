// Copy-number callers and tumour-fraction estimators, with the settings used in the cfCNV paper.
// Every output file is prefixed with "<sample_id>.<tool>.<config>" so all results can be collected together.

// Autosome names are read from the .fai at run time (1..22 or chr1..chr22) by this bash snippet:
//   CHRS = comma-separated autosomes, CHRSX = autosomes + X, STYLE = UCSC or NCBI
def CHR_SNIPPET(fai) {
    return """P=\$(head -n1 ${fai} | cut -f1 | grep -q '^chr' && echo chr || true)
    CHRS=\$(for i in \$(seq 1 22); do printf "%s%s," "\$P" "\$i"; done | sed 's/,\$//')
    CHRSX="\$CHRS,\${P}X"; STYLE=\$([ -n "\$P" ] && echo UCSC || echo NCBI)"""
}

process QDNASEQ {
    tag "${meta.id} ${bin_kb}kb"
    label 'cfcnv_r'

    input:
    tuple val(meta), path(bam), path(bai), val(bin_kb)

    output:
    tuple val(meta), val("${bin_kb}kb"), path("${meta.id}.QDNAseq.${bin_kb}kb.calls.vcf"), emit: calls
    tuple val(meta), val("${bin_kb}kb"), path("${meta.id}.QDNAseq.${bin_kb}kb.qdnaseq_seg.rds"), emit: rds
    path "${meta.id}.QDNAseq.${bin_kb}kb.segments.igv", emit: segments

    script:
    def g = params.genome == 'GRCh37' ? 'hg19' : 'hg38'
    """
    run_qdnaseq.R ${meta.id} ${bam} ${bin_kb} ${g} ${meta.id}.QDNAseq.${bin_kb}kb
    """
    stub:
    """
    printf '##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n1\\t10000001\\t.\\tN\\t<DEL>\\t.\\tPASS\\tSVTYPE=DEL;END=11000000\\n' > ${meta.id}.QDNAseq.${bin_kb}kb.calls.vcf
    touch ${meta.id}.QDNAseq.${bin_kb}kb.qdnaseq_seg.rds ${meta.id}.QDNAseq.${bin_kb}kb.segments.igv
    """
}

process ACE {
    tag "${meta.id} ${cfg}"
    label 'cfcnv_r'
    label 'small'

    input:
    tuple val(meta), val(cfg), path(rds)

    output:
    tuple val(meta), val(cfg), path("${meta.id}.ACE.${cfg}.ace.tsv"), emit: tf

    script:
    """
    run_ace.R ${meta.id} ${rds} ${meta.id}.ACE.${cfg}.ace.tsv
    """
    stub:
    """
    printf 'sample_id\\tace_cellularity_ploidy2\\tn_minima\\tace_cellularity_free\\tace_ploidy_free\\n${meta.id}\\t0.05\\t3\\t0.05\\t2\\n' > ${meta.id}.ACE.${cfg}.ace.tsv
    """
}

process ICHORCNA {
    tag "${meta.id} ${mode}"
    label 'cfcnv_r'
    label 'small'

    input:
    tuple val(meta), path(bam), path(bai), val(mode), val(cfg), val(tagc)
    path fai

    output:
    tuple val(meta), val(cfg), path("${meta.id}.ichorCNA.${tagc}.seg.txt"), emit: calls
    tuple val(meta), val(cfg), path("${meta.id}.ichorCNA.${tagc}.params.txt"), emit: tf
    path "${meta.id}.ichorCNA.${tagc}/*", emit: plots, optional: true

    script:
    def W    = params.ichor_window_kb as int
    def hg   = params.genome == 'GRCh37' ? 'hg19' : 'hg38'
    def pw   = W == 1000 ? '1Mb' : "${W}kb"
    def cent = params.genome == 'GRCh37' ? 'GRCh37.p13_centromere_UCSC-gapTable.txt' : 'GRCh38.GCA_000001405.2_centromere_acen.txt'
    def pon  = params.genome == 'GRCh37' ? "HD_ULP_PoN_${pw}_median_normAutosome_mapScoreFiltered_median.rds" : "HD_ULP_PoN_hg38_${pw}_median_normAutosome_median.rds"
    def pa   = mode == 'lowTF' ?
        "--ploidy 'c(2)' --normal 'c(0.95,0.99,0.995,0.999)' --maxCN 3 --estimateScPrevalence False --scStates 'c()'" :
        "--ploidy 'c(2,3)' --normal 'c(0.5,0.6,0.7,0.8,0.9,0.95,0.99)' --maxCN 5 --estimateScPrevalence True --scStates 'c(1,3)'"
    """
    ${CHR_SNIPPET(fai)}
    E=\${ICHOR_DIR:-/opt/ichorCNA}/inst/extdata
    readCounter --window ${W * 1000} --quality 20 --chromosome "\$CHRSX" ${bam} > ${meta.id}.wig
    Rscript \${ICHOR_DIR:-/opt/ichorCNA}/scripts/runIchorCNA.R --id ${meta.id} --WIG ${meta.id}.wig \\
        --gcWig \$E/gc_${hg}_${W}kb.wig --mapWig \$E/map_${hg}_${W}kb.wig --centromere \$E/${cent} \\
        --normalPanel \$E/${pon} ${pa} \\
        --includeHOMD False --chrs 'c(1:22, "X")' --chrTrain 'c(1:22)' \\
        --estimateNormal True --estimatePloidy True --txnE 0.9999 --txnStrength 10000 \\
        --genomeBuild ${hg} --genomeStyle \$STYLE --plotFileType png --outDir ichor_out
    cp ichor_out/${meta.id}.seg.txt    ${meta.id}.ichorCNA.${tagc}.seg.txt
    cp ichor_out/${meta.id}.params.txt ${meta.id}.ichorCNA.${tagc}.params.txt
    mkdir -p ${meta.id}.ichorCNA.${tagc} && cp -r ichor_out/${meta.id}/* ${meta.id}.ichorCNA.${tagc}/ 2>/dev/null || true
    """
    stub:
    """
    printf 'ID\\tchrom\\tstart\\tend\\tCorrected_Copy_Number\\n${meta.id}\\t1\\t10000001\\t12000000\\t1\\n' > ${meta.id}.ichorCNA.${tagc}.seg.txt
    printf 'Tumor Fraction:\\t0.012\\nPloidy:\\t2\\n' > ${meta.id}.ichorCNA.${tagc}.params.txt
    """
}

process CNVKIT_REFERENCE {
    // pooled-normal reference, built once when --cnvkit_normals is given
    tag "pooled ${params.cnvkit_bin}"
    label 'cfcnv_cnvkit'
    label 'medium'

    input:
    path normals
    tuple path(fasta), path(fai)

    output:
    path "cnvkit_reference.cnn", emit: cnn

    script:
    def n = normals.findAll { it.name.endsWith('.bam') }.join(' ')
    """
    cnvkit.py batch -m wgs -n ${n} -f ${fasta} --target-avg-size ${params.cnvkit_bin} \\
        --output-reference cnvkit_reference.cnn -d ref_build -p ${task.cpus}
    """
    stub:
    """
    touch cnvkit_reference.cnn
    """
}

process CNVKIT {
    tag "${meta.id}"
    label 'cfcnv_cnvkit'
    label 'medium'

    input:
    tuple val(meta), path(bam), path(bai)
    path cnn          // pooled reference, or NO_FILE -> flat reference (paper's real-data setting)
    tuple path(fasta), path(fai)
    val cfg

    output:
    tuple val(meta), val(cfg), path("${meta.id}.CNVkit.${cfg}.call.cns"), emit: calls
    path "${meta.id}.CNVkit.${cfg}.cnr", emit: cnr

    script:
    def stem = bam.baseName
    def ref  = cnn.name != 'NO_FILE' ? "-r ${cnn}" : "-n -f ${fasta} --target-avg-size ${params.cnvkit_bin}"
    """
    cnvkit.py batch ${bam} -m wgs ${ref} -d out -p ${task.cpus}
    cp out/${stem}.call.cns ${meta.id}.CNVkit.${cfg}.call.cns
    cp out/${stem}.cnr      ${meta.id}.CNVkit.${cfg}.cnr
    """
    stub:
    """
    printf 'chromosome\\tstart\\tend\\tgene\\tlog2\\tcn\\n1\\t10000000\\t11000000\\t-\\t-0.9\\t1\\n' > ${meta.id}.CNVkit.${cfg}.call.cns
    touch ${meta.id}.CNVkit.${cfg}.cnr
    """
}

process MOSDEPTH_BINS {
    tag "${meta.id}"
    label 'cfcnv_tools'
    label 'small'
    cpus 3

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.regions.bed.gz"), emit: bed

    script:
    """
    mosdepth -t ${task.cpus} -n -x -b ${params.spectre_bin} -Q 20 ${meta.id} ${bam}
    """
    stub:
    """
    echo | gzip > ${meta.id}.regions.bed.gz
    """
}

process SPECTRE {
    tag "${meta.id}"
    label 'cfcnv_spectre'
    label 'small'
    cpus 3

    input:
    tuple val(meta), path(cov)
    tuple path(fasta), path(fai)
    path blacklist
    val cfg

    output:
    tuple val(meta), val(cfg), path("${meta.id}.Spectre.${cfg}_cnv.bed.gz"), emit: calls

    script:
    def B = params.spectre_bin as int
    def minlen = Math.max(B, 50000)
    def bl = blacklist.name != 'NO_FILE' ? "--blacklist ${blacklist}" : ''
    """
    ${CHR_SNIPPET(fai)}
    spectre CNVCaller --coverage ${cov} --sample-id ${meta.id} --output-dir calls --reference ${fasta} \\
        --min-cnv-len ${minlen} --threads ${task.cpus} --cancer --only-chr "\$CHRS" ${bl}
    cp calls/${meta.id}_cnv.bed.gz ${meta.id}.Spectre.${cfg}_cnv.bed.gz
    """
    stub:
    """
    printf '1\\t10000000\\t11000000\\tDEL\\n' | gzip > ${meta.id}.Spectre.${cfg}_cnv.bed.gz
    """
}

process BICSEQ2_RESOURCES {
    // ONE-OFF, slow (~1-2 days): per-chromosome FASTA + GEM mappability for the read length
    tag "${params.bicseq2_readlen}-mer"
    label 'cfcnv_bicseq2'
    storeDir "${params.outdir}/resources/bicseq2"

    input:
    tuple path(fasta), path(fai)

    output:
    path "bicseq2_${params.bicseq2_readlen}mer", emit: dir

    script:
    def L = params.bicseq2_readlen
    """
    bicseq2_resources.sh ${fasta} ${L} bicseq2_${L}mer ${task.cpus}
    """
    stub:
    """
    mkdir -p bicseq2_${params.bicseq2_readlen}mer/reference
    """
}

process BICSEQ2_UNIQUE {
    tag "${meta.id}"
    label 'cfcnv_bicseq2'
    label 'medium'

    input:
    tuple val(meta), path(bam), path(bai)
    path fai

    output:
    tuple val(meta), path("${meta.id}_unique"), emit: unique

    script:
    """
    ${CHR_SNIPPET(fai)}
    mkdir -p ${meta.id}_unique
    H=\${BICSEQ2_GETUNIQUE:-/samtools-0.1.7a_getUnique-0.1.3/misc/samtools.pl}
    i=0
    for C in \$(echo "\$CHRS" | tr ',' ' '); do
      ( set -o pipefail; samtools view ${bam} \$C | perl \$H unique - | cut -f4 > ${meta.id}_unique/${meta.id}.\$C.seq ) &
      i=\$((i+1)); (( i % ${task.cpus} == 0 )) && wait
    done; wait
    ( set +o pipefail; samtools view -f 2 -F 3852 ${bam} | awk '\$9>0 && \$9<2000 {print \$9; n++; if(n>=100000) exit}' ) \\
      | sort -n | awk '{a[NR]=\$1} END{ if(NR<100) {print 167; exit}; print (NR%2)? a[(NR+1)/2] : int((a[NR/2]+a[NR/2+1])/2+0.5)}' \\
      > ${meta.id}_unique/fragment_size.txt
    """
    stub:
    """
    mkdir -p ${meta.id}_unique && echo 167 > ${meta.id}_unique/fragment_size.txt
    """
}

process BICSEQ2_NORMSEG {
    tag "${meta.id}"
    label 'cfcnv_bicseq2'

    input:
    tuple val(meta), path(unique)
    path resources
    path fai
    val cfg

    output:
    tuple val(meta), val(cfg), path("${meta.id}.BIC-seq2.${cfg}.cnv"), emit: calls

    script:
    def B = params.bicseq2_bin as int
    def L = params.bicseq2_readlen
    def mer = params.bicseq2_map_prefix ?: "bicseq2_${L}mer"
    """
    ${CHR_SNIPPET(fai)}
    R=\$(readlink -f ${resources}); U=\$(readlink -f ${unique}); mkdir -p norm seg tmp segtmp
    printf "chromName\\tfaFile\\tMapFile\\treadPosFile\\tbinFileNorm\\n" > norm.cfg
    printf "chromName\\tbinFileNorm\\n" > seg.cfg
    for C in \$(echo "\$CHRS" | tr ',' ' '); do
      printf "%s\\t%s\\t%s\\t%s\\t%s\\n" \$C \$R/reference/\$C.fa \$R/mappability/${mer}/${mer}.\$C.txt \$U/${meta.id}.\$C.seq \$PWD/norm/\$C.norm.bin >> norm.cfg
      printf "%s\\t%s\\n" \$C \$PWD/norm/\$C.norm.bin >> seg.cfg
    done
    perl \${BICSEQ2_NORM:-/NBICseq-norm_v0.2.4/NBICseq-norm.pl} --tmp=tmp -l ${L} -s \$(cat \$U/fragment_size.txt) \\
        -b ${B} --fig ${meta.id}.GC.pdf norm.cfg norm_parameters.txt
    perl \${BICSEQ2_SEG:-/NBICseq-seg_v0.7.2/NBICseq-seg.pl} --lambda=${params.bicseq2_lambda} --tmp=segtmp \\
        --fig=${meta.id}.seg.png seg.cfg ${meta.id}.BIC-seq2.${cfg}.cnv
    """
    stub:
    """
    printf 'chrom\\tstart\\tend\\tbinNum\\tobserved\\texpected\\tlog2.copyRatio\\n1\\t10000000\\t11000000\\t100\\t50\\t100\\t-1\\n' > ${meta.id}.BIC-seq2.${cfg}.cnv
    """
}

process FRAGLE {
    tag "${meta.id}"
    label 'cfcnv_fragle'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), val('raw WGS'), path("${meta.id}.Fragle.csv"), emit: tf

    script:
    def g = params.genome == 'GRCh37' ? 'hg19' : 'hg38'
    """
    mkdir -p in && ln -s \$(readlink -f ${bam}) in/${meta.id}.bam && ln -s \$(readlink -f ${bai}) in/${meta.id}.bam.bai
    cd \${FRAGLE_DIR:-/opt/FRAGLE}
    python main.py --input \$OLDPWD/in --output \$OLDPWD/fragle_out --mode R --genome_build ${g} --cpu ${task.cpus}
    cd -
    cp fragle_out/Fragle.csv ${meta.id}.Fragle.csv
    """
    stub:
    """
    printf 'Sample_ID,ctDNA_Burden\\n${meta.id},0.2\\n' > ${meta.id}.Fragle.csv
    """
}
