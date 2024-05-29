process GCPARAGON {
    tag "$meta.sample_id"
    label 'process_high'

    //container = 'ghcguzman/samtools1.20'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*/*.GCtagged.bam"), path("*/*.GCtagged.bam.bai"), emit: bam_bai

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    gcparagon \\
    --bam $bam \\
    --output-bam \\
    -rgb ${params.rgb} \\
    --threads $task.cpus
    """
}