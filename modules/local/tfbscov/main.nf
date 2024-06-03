process CREATE_TFBSCOV_MATRIX {
    tag "$meta.sample_id"
    label 'process_medium'

    //container = 'ghcguzman/samtools1.20'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.1Mb_bins.tsv"), emit: matrix

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"
    """
    echo "hi"
    """
}