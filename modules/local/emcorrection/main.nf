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
    -t $genome \\
    -b $blacklist \\
    -obam ${bam.baseName}.EMtagged.bam \\
    -n $task.cpus \\
    -cw ${bam.baseName}.EMtagged.em.weights.tsv \\
    -fcw ${bam.baseName}.EMtagged.fl.weights.tsv \\
    -cm ${params.em_norm} \\
    -l ${bam.baseName}.EMtagged.log
    """
}