# GENECE WISHBONE Pipeline

Take a TSV of samples, apply fragment-based GC correction, then generate GC-corrected COVERAGE and FEMS matrixes.

## Pipeline Steps

1. Compute genome-wide fragment-based GC bias for each sample and correct (GCparagon).
2. Compute GC-corrected COVERAGE matrix (custom).
3. Compute GC-corrected FEMS matrix (custom).

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

**NOTE:** The `--gc_correction` **AND** `--em_correction` arguments are required if you want to run end motif correction on your bam files.

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

1. Create docker containers for each step to ensure reproducibility.

2. Implement TFBSCov.

## CHANGELOG

`May 29, 2024` : Created first version of the `WISHBONE` pipeline.

`June 02, 2024` : Implemented `WISHBONE` pipeline `v2`. Changed output of COVERAGE matrix to reflect GENECE pipeline input. Added preliminary end-motif correction. Added `--gc_correction` and `--em_correction` parameters.