process EMCORRECTION {
    tag "$meta.sample_id"
    label 'process_highest'

    container = 'ghcguzman/emcorrection.amd64'

    input:
    tuple val(meta), path(bam), path(bai)
    path(genome)
    path(blacklist)

    output:
    tuple val(meta), path("*.EMtagged.bam"), emit: bam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    emcorrection.py \\
    -i $bam \\
    -r $genome \\
    -e $blacklist \\
    -o ${bam.baseName}.EMtagged.bam \\
    --threads $task.cpus \\
    --sort_bam
    """
}