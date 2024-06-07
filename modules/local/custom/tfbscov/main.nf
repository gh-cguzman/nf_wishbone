process CREATE_TFBSCOV_MATRIX {
    tag "$meta.sample_id"
    label 'process_low'

    //container = 'ghcguzman/samtools1.20'

    input:
    tuple val(meta), path(bam), path(bai)
    path(motif_beds)

    output:
    tuple val(meta), path("*.mat.smoothed.norm.tsv")     , emit: smooth_norm_mat
    tuple val(meta), path("*.features.smoothed.norm.tsv"), emit: smooth_norm_features
    tuple val(meta), path("plots/*coverage_plot.png")   , emit: cov_png
    tuple val(meta), path("plots/*fft_plot.png")        , emit: fft_png

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"
    """
    mkdir -p plots/

    bam2tfbscov.py \\
    -B $bam \\
    -b $motif_beds \\
    -o ${bam.baseName}.mat.smoothed.norm.tsv \\
    -s ${bam.baseName}

    tfbsmat2features.py \\
    --input ${bam.baseName}.mat.smoothed.norm.tsv \\
    --output ${bam.baseName}.features.smoothed.norm.tsv \\
    --plot_dir plots/
    """
}