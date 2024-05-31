process CREATE_COVERAGE_MATRIX {
    tag "$meta.sample_id"
    label 'process_medium'

    //container = 'ghcguzman/samtools1.20'

    input:
    tuple val(meta), path(bam), path(bai)
    path(regions)
    path(blacklist)

    output:
    tuple val(meta), path("*.1Mb_bins.tsv"), emit: matrix

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"
    """
    bam2cov.py \\
    --bam $bam \\
    --regions $regions \\
    --blacklist $blacklist \\
    --outcounts ${bam.baseName}.1Mb_bins.tsv.tmp \\
    --outlog ${bam.baseName}.log \\
    --threads $task.cpus \\
    --sample ${bam.baseName}

    convert.py \\
    --cov ${bam.baseName}.1Mb_bins.tsv.tmp \\
    --sample_id ${meta.sample_id} \\
    --output ${bam.baseName}.1Mb_bins.tsv
    """
}