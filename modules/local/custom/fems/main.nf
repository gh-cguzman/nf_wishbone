process CREATE_FEMS_MATRIX {
    tag "$meta.sample_id"
    label 'process_medium'

    //container = 'ghcguzman/samtools1.20'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.FEMS.df.tsv"), emit: matrix

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"
    """
    bam2fems.py \\
    --bam $bam \\
    --out_em ${bam.baseName}.EM.df.tsv \\
    --out_fems ${bam.baseName}.FEMS.df.tsv \\
    --out_qc ${bam.baseName}.QC.df.tsv
    """
}