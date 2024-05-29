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

### Starting with **NON** GC-corrected BAM files
```
nextflow run gh-cguzman/nf_wishbone -r main -profile macbook,docker --input input.tsv
```

### Starting with GC-corrected BAM files
```
nextflow run gh-cguzman/nf_wishbone -r main -profile macbook,docker --input input.tsv --gc_corrected
```

**NOTE:** The `--gc_corrected` argument is required if you want to avoid running `gcparagon` on the bam files and proceed directly to feature extraction.

### Running the pipeline on the `.88` server
```
nextflow run gh-cguzman/nf_tfbscov -r main -profile server88 --input input.tsv
```

**NOTE:** This command assumes that you are working in an environment that has all the necessary packages installed. Docker is not available on the `.88`.

### Running the pipeline on the `GCG` server
```
nextflow run gh-cguzman/nf_tfbscov -r main -profile robin14 --input input.tsv
```

**NOTE:** This command assumes that you are working in an environment that has all the necessary packages installed. Docker is not available on GCG's `robin14`.

## TODO

## CHANGELOG

`May 29, 2024` : Created first version of the `WISHBONE` pipeline.