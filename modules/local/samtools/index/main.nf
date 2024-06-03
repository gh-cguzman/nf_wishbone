process SAMTOOLS_INDEX {
    tag "$meta.sample_id"
    label 'process_single'

    conda (params.enable_conda ? 'bioconda::samtools=1.20' : null)
    container = 'ghcguzman/samtools1.20'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam", includeInputs: true), path("*.bai"), emit: bam_bai

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools index \\
    $bam
    """
}