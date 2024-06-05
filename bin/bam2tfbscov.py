#!/usr/bin/env python
import argparse
import pysam
import pybedtools
import numpy as np
import pandas as pd
import logging
from scipy.signal import savgol_filter
from multiprocessing import Pool, cpu_count
from functools import partial
import os

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_arguments():
    parser = argparse.ArgumentParser(description="Calculate GC tag sums across genomic positions for multiple BED files")
    parser.add_argument("--bam", required=True, help="Input GC corrected BAM file")
    parser.add_argument("--beds", required=True, nargs='+', help="Input BED files of genomic positions")
    parser.add_argument("--min_fl", type=int, default=110, help="Minimum fragment length")
    parser.add_argument("--max_fl", type=int, default=230, help="Maximum fragment length")
    parser.add_argument("--output", required=True, help="Output file for the normalized and smoothed matrix")
    parser.add_argument("--output_raw", required=True, help="Output file for the raw matrix")
    parser.add_argument("--threads", type=int, default=cpu_count(), help="Number of threads to use for parallel processing")
    parser.add_argument("--range", type=int, default=1000, help="Range around center to consider for GC tag sums (+/-)")
    parser.add_argument("--sample_id", required=True, help="Sample ID for the output file")
    return parser.parse_args()

def get_gc_tag_sums(bamfile, chrom, center, min_fl, max_fl, window_range):
    gc_sums = np.zeros(window_range * 2 + 1)
    start = center - window_range
    end = center + window_range

    if chrom not in [f'chr{i}' for i in range(1, 23)]:
        return gc_sums

    for read in bamfile.fetch(chrom, start, end):
        if not read.is_proper_pair or read.is_unmapped or read.mapping_quality == 0 or read.is_supplementary:
            continue
        if min_fl <= read.template_length <= max_fl:
            gc_pos = read.reference_start - start
            if 0 <= gc_pos < len(gc_sums):
                gc_sums[gc_pos] += read.get_tag("GC")

    return gc_sums

def calculate_normalization_factors(bamfile, chrom, center, min_fl, max_fl):
    upstream_start = center - 1500
    upstream_end = center
    downstream_start = center
    downstream_end = center + 1500

    upstream_gc_sum = 0
    downstream_gc_sum = 0

    if chrom not in [f'chr{i}' for i in range(1, 23)]:
        return 1

    for read in bamfile.fetch(chrom, upstream_start, upstream_end):
        if not read.is_proper_pair or read.is_unmapped or read.mapping_quality == 0 or read.is_supplementary:
            continue
        if min_fl <= read.template_length <= max_fl:
            upstream_gc_sum += read.get_tag("GC")

    for read in bamfile.fetch(chrom, downstream_start, downstream_end):
        if not read.is_proper_pair or read.is_unmapped or read.mapping_quality == 0 or read.is_supplementary:
            continue
        if min_fl <= read.template_length <= max_fl:
            downstream_gc_sum += read.get_tag("GC")

    normalization_factor = (upstream_gc_sum + downstream_gc_sum) / 3000
    if normalization_factor == 0:
        normalization_factor = 1  # Avoid division by zero

    return normalization_factor

def process_bed_file(bed_file, args):
    bed = pybedtools.BedTool(bed_file)
    base_name = os.path.basename(bed_file).split('.')[0]
    min_fl = args.min_fl
    max_fl = args.max_fl
    window_range = args.range

    bamfile = pysam.AlignmentFile(args.bam, "rb")

    logging.info(f"Processing intervals for {bed_file} ...")
    intervals = [(interval.chrom, (interval.start + interval.end) // 2) for interval in bed]

    normalized_gc_sums_list = np.zeros((len(intervals), window_range * 2 + 1))
    raw_gc_sums_list = np.zeros((len(intervals), window_range * 2 + 1))

    for i, (chrom, center) in enumerate(intervals):
        gc_sums = get_gc_tag_sums(bamfile, chrom, center, min_fl, max_fl, window_range)
        normalization_factor = calculate_normalization_factors(bamfile, chrom, center, min_fl, max_fl)

        normalized_gc_sums = gc_sums / normalization_factor
        normalized_gc_sums_list[i] = normalized_gc_sums
        raw_gc_sums_list[i] = gc_sums

    bamfile.close()

    # Remove bins with extremely high coverage (10 standard deviations above the mean)
    logging.info(f"Removing outlier bins for {bed_file} ...")
    mean_coverage = np.mean(normalized_gc_sums_list, axis=0)
    std_coverage = np.std(normalized_gc_sums_list, axis=0)
    filtered_gc_sums = normalized_gc_sums_list[np.all(normalized_gc_sums_list <= mean_coverage + 10 * std_coverage, axis=1)]

    mean_gc_sums = np.mean(filtered_gc_sums, axis=0)

    # Smooth coverage profiles with window length 165 and polynomial order 3
    smoothed_gc_sums = savgol_filter(mean_gc_sums, window_length=165, polyorder=3)

    # Normalize coverage profile to mean coverage of 1 across +/-1000bp window
    normalized_gc_sums = smoothed_gc_sums / np.mean(smoothed_gc_sums)

    # Create the single-row output for normalized and smoothed data
    output_row = [args.sample_id, base_name] + normalized_gc_sums.tolist()

    # Create the single-row output for raw data
    raw_output_row = [args.sample_id, base_name] + np.mean(raw_gc_sums_list, axis=0).tolist()

    return output_row, raw_output_row

def process_bed_files_in_parallel(args, bed_files):
    process_func = partial(process_bed_file, args=args)
    chunksize = max(1, len(bed_files) // (2 * args.threads))  # Dynamic chunksize calculation
    with Pool(args.threads) as pool:
        results = list(pool.imap(process_func, bed_files, chunksize=chunksize))
    return results

def main():
    args = parse_arguments()

    window_range = args.range
    results = []
    raw_results = []

    logging.info("Processing BED files in parallel...")
    parallel_results = process_bed_files_in_parallel(args, args.beds)

    logging.info("Staging results ...")
    for result_pair in parallel_results:
        results.append(result_pair[0])
        raw_results.append(result_pair[1])

    # Save the results to output files
    positions = np.arange(-window_range, window_range + 1)
    output_columns = ["Sample_ID", "BED_File"] + [str(i) for i in positions]
    output_df = pd.DataFrame(results, columns=output_columns)
    output_df.to_csv(args.output, sep='\t', index=False)

    raw_output_df = pd.DataFrame(raw_results, columns=output_columns)
    raw_output_df.to_csv(args.output_raw, sep='\t', index=False)

    logging.info("Script completed successfully")

if __name__ == "__main__":
    main()
