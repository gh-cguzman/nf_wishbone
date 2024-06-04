import argparse
import pysam
import pybedtools
import numpy as np
import pandas as pd
import logging
from scipy.signal import savgol_filter
from concurrent.futures import ProcessPoolExecutor, as_completed

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
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use for parallel processing")
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

def process_interval(args):
    bamfile_path, chrom, center, min_fl, max_fl, window_range = args
    bamfile = pysam.AlignmentFile(bamfile_path, "rb")
    gc_sums = get_gc_tag_sums(bamfile, chrom, center, min_fl, max_fl, window_range)
    bamfile.close()
    return gc_sums

def calculate_normalization_factors(bamfile, chrom, center, min_fl, max_fl):
    upstream_start = center - 3000
    upstream_end = center - 1500
    downstream_start = center + 1500
    downstream_end = center + 3000

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

    normalization_factor = (upstream_gc_sum + downstream_gc_sum) / (3000 * 2)
    if normalization_factor == 0:
        normalization_factor = 1  # Avoid division by zero

    return normalization_factor

def process_bed_file(bed_file, args):
    bed = pybedtools.BedTool(bed_file)
    base_name = bed_file.split('/')[-1].split('.')[0]
    bamfile_path = args.bam
    min_fl = args.min_fl
    max_fl = args.max_fl
    window_range = args.range

    logging.info(f"Processing intervals for {bed_file} ...")
    intervals = [(bamfile_path, interval.chrom, (interval.start + interval.end) // 2, min_fl, max_fl, window_range) for interval in bed]

    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        all_gc_sums = list(executor.map(process_interval, intervals))

    normalized_gc_sums_list = []
    raw_gc_sums_list = []

    logging.info(f"Calculating normalization factors for {bed_file} ...")
    for interval, gc_sums in zip(intervals, all_gc_sums):
        chrom = interval[1]
        center = interval[2]
        bamfile = pysam.AlignmentFile(bamfile_path, "rb")
        normalization_factor = calculate_normalization_factors(bamfile, chrom, center, min_fl, max_fl)
        bamfile.close()

        normalized_gc_sums = gc_sums / normalization_factor
        normalized_gc_sums_list.append(normalized_gc_sums)
        raw_gc_sums_list.append(gc_sums)

    # Remove bins with extremely high coverage (10 standard deviations above the mean)
    logging.info(f"Removing outlier bins for {bed_file} ...")
    normalized_gc_sums_array = np.array(normalized_gc_sums_list)
    mean_coverage = np.mean(normalized_gc_sums_array, axis=0)
    std_coverage = np.std(normalized_gc_sums_array, axis=0)
    filtered_gc_sums = np.array([gc_sums for gc_sums in normalized_gc_sums_array if not np.any(gc_sums > mean_coverage + 10 * std_coverage)])

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

def main():
    args = parse_arguments()

    logging.info(f"Loading BAM {args.bam.split('/')[-1]} file ...")
    bamfile = pysam.AlignmentFile(args.bam, "rb")
    bamfile.close()

    window_range = args.range
    results = []
    raw_results = []

    with ProcessPoolExecutor(max_workers=len(args.beds)) as executor:
        futures = {executor.submit(process_bed_file, bed_file, args): bed_file for bed_file in args.beds}
        for future in as_completed(futures):
            bed_file = futures[future]
            try:
                result, raw_result = future.result()
                results.append(result)
                raw_results.append(raw_result)
            except Exception as e:
                logging.error(f"Error processing {bed_file}: {e}")

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
