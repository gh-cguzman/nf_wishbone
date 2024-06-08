#!/usr/bin/env python3
import argparse
import logging
import pysam
import numpy as np
from scipy.signal import savgol_filter
from concurrent.futures import ProcessPoolExecutor, as_completed

def setup_logging():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_bed_file(bed_file):
    regions = []
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('track') or line.startswith('browser'):
                continue
            chrom, start, end = line.strip().split()[:3]
            regions.append((chrom, int(start), int(end)))
    return regions

def extend_bed_regions(bed, window=5000):
    extended_bed = [(chrom, max(0, start - window), end + window) for chrom, start, end in bed]
    return extended_bed

def filter_chr(bed):
    chr_filter = {f'chr{i}' for i in range(1, 23)}
    return [(chrom, start, end) for chrom, start, end in bed if chrom in chr_filter]

def human_sort_bed(bed):
    chrom_order = {f'chr{i}': i for i in range(1, 23)}
    return sorted(bed, key=lambda x: (chrom_order[x[0]], x[1], x[2]))

def process_bam(bam_file, bed_chunk, min_fragment_size, max_fragment_size):
    logging.info(f"Processing BAM file: {bam_file}")
    regions_gc_content = []

    bamfile = pysam.AlignmentFile(bam_file, "rb")

    for chrom, start, end in bed_chunk:
        array = np.zeros(end - start + 1)  # Initialize the array to zeros

        for read in bamfile.fetch(chrom, start, end):
            if not read.is_proper_pair or read.is_duplicate or read.is_supplementary:
                continue
            if min_fragment_size <= read.template_length <= max_fragment_size:
                midpoint = read.reference_start + read.template_length // 2
                if read.has_tag('GC'):
                    gc_value = read.get_tag('GC')
                    pos_in_array = midpoint - start
                    if 0 <= pos_in_array < len(array):
                        array[pos_in_array] += gc_value

        regions_gc_content.append((chrom, start, end, array))

    bamfile.close()

    return regions_gc_content

def filter_high_coverage_bins(regions_gc_content, threshold=10):
    logging.info(f"Filtering out regions with coverage above {threshold} standard deviations")
    all_values = np.concatenate([array for _, _, _, array in regions_gc_content])
    mean = np.mean(all_values)
    std = np.std(all_values)
    cutoff = mean + threshold * std

    filtered_regions = []
    for chrom, start, end, array in regions_gc_content:
        array[array > cutoff] = 0
        filtered_regions.append((chrom, start, end, array))

    return filtered_regions

def average_coverage_by_source(filtered_regions):
    logging.info("Computing average coverage for each source file")
    all_arrays = np.vstack([array for _, _, _, array in filtered_regions])
    avg_coverage = np.mean(all_arrays, axis=0)
    return avg_coverage

def smooth_and_normalize(average):
    logging.info("Smoothing and normalizing coverage")
    mean_coverage = np.mean(average)
    if (mean_coverage > 0) and (average.size >= 6001):
        normalized_coverage = average / mean_coverage
        coverage_region = normalized_coverage[4000:6001]
        filtered_array = savgol_filter(coverage_region, window_length=165, polyorder=3)
    else:
        filtered_array = np.zeros(2001)  # In case of invalid data
    return filtered_array

def write_averaged_coverage_to_tsv(output_file, sample_id, bed_file, avg_coverage):
    logging.info(f"Writing averaged coverage to TSV file: {output_file}")
    with open(output_file, 'a') as f:
        line = [sample_id, bed_file.split('.')[0]] + list(map(str, avg_coverage))
        f.write("\t".join(line) + "\n")

def process_bed_file_in_parallel(bed_file, bam_file, min_fragment_size, max_fragment_size, sample_id, chunk_size):
    bed = parse_bed_file(bed_file)
    extended_bed = extend_bed_regions(bed, window=5000)
    filtered_bed = filter_chr(extended_bed)
    sorted_bed = human_sort_bed(filtered_bed)

    all_processed_regions = []
    for i in range(0, len(sorted_bed), chunk_size):
        bed_chunk = sorted_bed[i:i + chunk_size]
        processed_regions = process_bam(bam_file, bed_chunk, min_fragment_size, max_fragment_size)
        all_processed_regions.extend(processed_regions)

    filtered_regions = filter_high_coverage_bins(all_processed_regions, threshold=10)

    averaged_coverage = average_coverage_by_source(filtered_regions)

    smoothed_normalized_averages = smooth_and_normalize(averaged_coverage)

    return bed_file, smoothed_normalized_averages

def main(args):
    setup_logging()

    bed_files = args.bed_files
    bam_file = args.bam_file
    output_file = args.output_file
    sample_id = args.sample_id
    min_fragment_size = args.min_fragment_size
    max_fragment_size = args.max_fragment_size
    num_cpus = args.num_cpus
    chunk_size = args.chunk_size

    with open(output_file, 'w') as f:
        header = ["Sample_ID", "BED_File"] + [str(i) for i in range(-1000, 1001)]
        f.write("\t".join(header) + "\n")

    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        futures = [executor.submit(process_bed_file_in_parallel, bed_file, bam_file, min_fragment_size, max_fragment_size, sample_id, chunk_size) for bed_file in bed_files]

        with open(output_file, 'a') as f:
            for future in as_completed(futures):
                bed_file, smoothed_normalized_averages = future.result()
                line = [sample_id, bed_file.split('.')[0]] + list(map(str, smoothed_normalized_averages))
                f.write("\t".join(line) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process BAM and BED files.")
    parser.add_argument("-b", "--bed-files", nargs='+', required=True, help="List of BED files.")
    parser.add_argument("-B", "--bam-file", required=True, help="Input BAM file.")
    parser.add_argument("-o", "--output-file", required=True, help="Output TSV file.")
    parser.add_argument("-s", "--sample-id", required=True, help="Sample ID.")
    parser.add_argument("-m", "--min-fragment-size", type=int, default=110, help="Minimum fragment size.")
    parser.add_argument("-M", "--max-fragment-size", type=int, default=230, help="Maximum fragment size.")
    parser.add_argument("-c", "--num-cpus", type=int, default=1, help="Number of CPUs to use for parallelization.")
    parser.add_argument("-cs", "--chunk-size", type=int, default=10000, help="Number of regions to process in a chunk.")

    args = parser.parse_args()
    main(args)
