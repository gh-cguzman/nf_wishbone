#!/usr/bin/env python3
import argparse
import logging
import pandas as pd
import pysam
import numpy as np
from scipy.signal import savgol_filter
import multiprocessing as mp

def setup_logging():
    """Setup logging configuration."""
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def extend_bed_regions(bed, window=5000):
    """Extend BED regions by ±window bp."""
    bed[1] = bed[1] - window
    bed[2] = bed[2] + window
    bed[1] = bed[1].clip(lower=0)  # Ensure start is not less than 0
    return bed

def filter_chr(bed):
    """Keep only regions in chr1 to chr22."""
    chr_filter = ['chr' + str(i) for i in range(1, 23)]
    bed_filtered = bed[bed[0].isin(chr_filter)].copy()
    return bed_filtered

def human_sort_bed(bed):
    """Sort BED file in human-readable order (chr1, chr2, ..., chr10, chr11, ..., chr22)."""
    bed.loc[:, 'chrom'] = pd.Categorical(bed[0], categories=[f'chr{i}' for i in range(1, 23)], ordered=True)
    sorted_bed = bed.sort_values(['chrom', 1, 2])
    return sorted_bed.drop(columns='chrom')

def process_bam(bam_file, bed, min_fragment_size, max_fragment_size):
    """Process BAM file and calculate the GC content in extended BED regions."""
    logging.info(f"Processing BAM file: {bam_file}")
    bed_df = bed.copy()
    bed_df['array'] = [np.zeros(end - start + 1) for start, end in zip(bed_df[1], bed_df[2])]
    bamfile = pysam.AlignmentFile(bam_file, "rb")

    for bed_region in bed_df.itertuples():
        chrom = bed_region[1]
        start = bed_region[2]
        end = bed_region[3]
        array = np.zeros(end - start + 1)  # Reset the array to zeros

        for read in bamfile.fetch(chrom, start, end):
            if not read.is_proper_pair or read.is_duplicate or read.is_supplementary:
                continue
            if min_fragment_size <= read.template_length <= max_fragment_size:
                midpoint = read.reference_start + read.template_length // 2
                gc_value = read.get_tag('GC')
                pos_in_array = midpoint - start
                if 0 <= pos_in_array < len(array):
                    array[pos_in_array] += gc_value

        bed_df.at[bed_region.Index, 'array'] = array  # Update the array in the DataFrame

    bamfile.close()

    return bed_df

def filter_high_coverage_bins(bed_df, threshold=10):
    """Filter out regions with extremely high coverage (above threshold standard deviations)."""
    logging.info(f"Filtering out regions with coverage above {threshold} standard deviations")
    all_values = np.concatenate(bed_df['array'].values)
    mean = np.mean(all_values)
    std = np.std(all_values)
    cutoff = mean + threshold * std

    bed_df['filtered_array'] = np.empty((len(bed_df),), dtype=object)
    for index, row in bed_df.iterrows():
        array = row['array']
        array[array > cutoff] = 0
        bed_df.at[index, 'filtered_array'] = array
    return bed_df

def average_coverage_by_source(bed_df):
    """Compute the average coverage for each source file."""
    logging.info("Computing average coverage for each source file")
    filtered_arrays = np.vstack(bed_df['filtered_array'].values)
    avg_coverage = np.mean(filtered_arrays, axis=0)
    return avg_coverage

def smooth_and_normalize(average):
    """Smooth and normalize arrays."""
    logging.info("Smoothing and normalizing coverage")
    # Normalize to a mean coverage of 1 across the -5000 to +5000 region
    mean_coverage = np.mean(average)
    if (mean_coverage > 0) and (average.size >= 6001):
        normalized_coverage = average / mean_coverage
        # Extract -1000 to +1000 region
        coverage_region = normalized_coverage[4000:6001]
        # Smooth the coverage region
        filtered_array = savgol_filter(coverage_region, window_length=165, polyorder=3)
    else:
        filtered_array = np.zeros(2001) # In case of invalid data
    return filtered_array

def save_averaged_coverage_to_tsv(averages, output_file, sample_id):
    """Save the averaged coverage data to a TSV file."""
    logging.info(f"Saving averaged coverage to TSV file: {output_file}")
    with open(output_file, 'w') as f:
        header = ["Sample_ID", "BED_File"] + [str(i) for i in range(-1000, 1001)]
        f.write("\t".join(header) + "\n")

        for bed_file, avg_coverage in averages.items():
            line = [sample_id, bed_file.split('.')[0]] + list(map(str, avg_coverage))
            f.write("\t".join(line) + "\n")

def worker_func(bed_file, bam_file, min_fragment_size, max_fragment_size, output_queue):
    """Worker function to process a single BED file."""
    logging.info(f"Processing BED file: {bed_file}")
    bed = pd.read_csv(bed_file, sep='\t', header=None, skiprows=1)
    extended_bed = extend_bed_regions(bed, window=5000)
    filtered_bed = filter_chr(extended_bed)
    sorted_bed = human_sort_bed(filtered_bed)

    processed_bed_df = process_bam(bam_file, sorted_bed, min_fragment_size, max_fragment_size)

    filtered_bed_df = filter_high_coverage_bins(processed_bed_df, threshold=10)

    averaged_coverage = average_coverage_by_source(filtered_bed_df)

    smoothed_normalized_averages = smooth_and_normalize(averaged_coverage)

    output_queue.put((bed_file, smoothed_normalized_averages))

def main(args):
    setup_logging()

    bed_files = args.bed_files
    bam_file = args.bam_file
    output_file = args.output_file
    sample_id = args.sample_id
    min_fragment_size = args.min_fragment_size
    max_fragment_size = args.max_fragment_size
    num_cpus = args.num_cpus

    output_queue = mp.Queue()
    print(f"Using {num_cpus} CPUs for parallelization.")

    processes = []

    for i, bed_file in enumerate(bed_files):
        if i >= num_cpus:
            break  # Don't create more processes than requested CPUs

        p = mp.Process(target=worker_func, args=(bed_file, bam_file, min_fragment_size, max_fragment_size, output_queue))
        processes.append(p)
        p.start()

    all_averages = {}
    for _ in range(len(bed_files)):
        bed_file, smoothed_normalized_averages = output_queue.get()
        all_averages[bed_file] = smoothed_normalized_averages

    for p in processes:
        p.join()

    save_averaged_coverage_to_tsv(all_averages, output_file, sample_id)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process BAM and BED files.")
    parser.add_argument("-b", "--bed-files", nargs='+', required=True, help="List of BED files.")
    parser.add_argument("-B", "--bam-file", required=True, help="Input BAM file.")
    parser.add_argument("-o", "--output-file", required=True, help="Output TSV file.")
    parser.add_argument("-s", "--sample-id", required=True, help="Sample ID.")
    parser.add_argument("-m", "--min-fragment-size", type=int, default=110, help="Minimum fragment size.")
    parser.add_argument("-M", "--max-fragment-size", type=int, default=230, help="Maximum fragment size.")
    parser.add_argument("-c", "--num-cpus", type=int, default=1, help="Number of CPUs to use for parallelization.")

    args = parser.parse_args()
    main(args)