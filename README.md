# GENECE WISHBONE Pipeline

Take a TSV of samples, apply fragment-based GC correction, then generate GC-corrected COVERAGE and FEMS matrixes.

## Pipeline Steps

1. Compute genome-wide fragment-based GC bias for each sample and correct (GCparagon).
2. Compute genome-wide fragment-based EM bias for each sample and correct (custom).
3. Compute bias-corrected COVERAGE features (custom).
4. Compute bias-corrected FEMS features (custom).
5. Compute bias-corrected TFBSCov features (custom).

## Updating the pipeline
```
nextflow pull gh-cguzman/nf_wishbone
```

## Running the pipeline

### Starting with **NON** corrected BAM files
```
nextflow run gh-cguzman/nf_wishbone -r main -profile macbook,docker --input input.tsv
```
OR
```
nextflow run gh-cguzman/nf_wishbone -r main -profile macbook,docker --bams '*.bam'
```

**NOTE:** The ''s are necessary for Nextflow to properly create channels of bam files.

### GC Correction
```
nextflow run gh-cguzman/nf_wishbone -r main -profile macbook,docker --input input.tsv --gc_correction
```

**NOTE:** The `--gc_correction` argument is required if you want to run `gcparagon` on your bam files.

### End Motif Correction
```
nextflow run gh-cguzman/nf_wishbone -r main -profile macbook,docker --input input.tsv --gc_correction --em_correction
```

**NOTE:** The `--gc_correction` **AND** `--em_correction` arguments are required if you want to run gc correction and end motif correction respectively on your bam files. Setting only one or the other will only carry out their respective corrections. By default the pipeline assumes that the input is **ALREADY** bias corrected and will do neither.

### Running the pipeline on the `.88` server
```
nextflow run gh-cguzman/nf_wishbone -r main -profile server88 --input input.tsv
```

**NOTE:** This command assumes that you are working in an environment that has all the necessary packages installed. Docker is not available on the `.88`.

### Running the pipeline on the `GCG` server
```
nextflow run gh-cguzman/nf_wishbone -r main -profile robin14 --input input.tsv
```

**NOTE:** This command assumes that you are working in an environment that has all the necessary packages installed.

## TODO

1. Create docker containers for each step to ensure reproducibility. There are issues here with Docker permissions on GCG server, and Singularity not liking `arm64` Docker images.

## CHANGELOG

`May 29, 2024` : Created first version of the `WISHBONE` pipeline.

`June 02, 2024` : Implemented `WISHBONE` pipeline `v2`. Changed output of COVERAGE matrix to reflect GENECE pipeline input. Added preliminary end-motif correction. Added `--gc_correction` and `--em_correction` parameters.

`June 04, 2024` : Added preliminary `TFBSCov` feature workflow.