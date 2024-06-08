#!/usr/bin/env python3
import argparse
import logging
import pysam
import numpy as np
from scipy.signal import savgol_filter
import multiprocessing as mp

def setup_logging():
    """Setup logging configuration."""
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_bed_file(bed_file):
    """Parse BED file into a list of tuples (chrom, start, end)."""
    regions = []
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('track') or line.startswith('browser'):
                continue
            chrom, start, end = line.strip().split()[:3]
            regions.append((chrom, int(start), int(end)))
    return regions

def extend_bed_regions(bed, window=5000):
    """Extend BED regions by ±window bp."""
    extended_bed = [(chrom, max(0, start - window), end + window) for chrom, start, end in bed]
    return extended_bed

def filter_chr(bed):
    """Keep only regions in chr1 to chr22."""
    chr_filter = {f'chr{i}' for i in range(1, 23)}
    return [(chrom, start, end) for chrom, start, end in bed if chrom in chr_filter]

def human_sort_bed(bed):
    """Sort BED file in human-readable order (chr1, chr2, ..., chr10, chr11, ..., chr22)."""
    chrom_order = {f'chr{i}': i for i in range(1, 23)}
    return sorted(bed, key=lambda x: (chrom_order[x[0]], x[1], x[2]))

def process_bam(bam_file, bed, min_fragment_size, max_fragment_size):
    """Process BAM file and calculate the GC content in extended BED regions."""
    logging.info(f"Processing BAM file: {bam_file}")
    regions_gc_content = []

    bamfile = pysam.AlignmentFile(bam_file, "rb")

    for chrom, start, end in bed:
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
    """Filter out regions with extremely high coverage (above threshold standard deviations)."""
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
    """Compute the average coverage for each source file."""
    logging.info("Computing average coverage for each source file")
    all_arrays = np.vstack([array for _, _, _, array in filtered_regions])
    avg_coverage = np.mean(all_arrays, axis=0)
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

def write_averaged_coverage_to_tsv(output_file, sample_id, bed_file, avg_coverage):
    """Append the averaged coverage data to a TSV file."""
    logging.info(f"Writing averaged coverage to TSV file: {output_file}")
    with open(output_file, 'a') as f:
        line = [sample_id, bed_file.split('.')[0]] + list(map(str, avg_coverage))
        f.write("\t".join(line) + "\n")

def worker_func(bed_file, bam_file, min_fragment_size, max_fragment_size, sample_id, queue, semaphore):
    """Worker function to process a single BED file and write the result to the output file."""
    with semaphore:
        logging.info(f"Processing BED file: {bed_file}")
        bed = parse_bed_file(bed_file)
        extended_bed = extend_bed_regions(bed, window=5000)
        filtered_bed = filter_chr(extended_bed)
        sorted_bed = human_sort_bed(filtered_bed)

        processed_regions = process_bam(bam_file, sorted_bed, min_fragment_size, max_fragment_size)

        filtered_regions = filter_high_coverage_bins(processed_regions, threshold=10)

        averaged_coverage = average_coverage_by_source(filtered_regions)

        smoothed_normalized_averages = smooth_and_normalize(averaged_coverage)

        queue.put((bed_file, smoothed_normalized_averages))

def main(args):
    setup_logging()

    bed_files = args.bed_files
    bam_file = args.bam_file
    output_file = args.output_file
    sample_id = args.sample_id
    min_fragment_size = args.min_fragment_size
    max_fragment_size = args.max_fragment_size
    num_cpus = args.num_cpus

    # Create the output file and write the header
    with open(output_file, 'w') as f:
        header = ["Sample_ID", "BED_File"] + [str(i) for i in range(-1000, 1001)]
        f.write("\t".join(header) + "\n")

    print(f"Using {num_cpus} CPUs for parallelization.")

    queue = mp.Queue()
    semaphore = mp.Semaphore(num_cpus)

    def start_worker(bed_file):
        p = mp.Process(target=worker_func, args=(bed_file, bam_file, min_fragment_size, max_fragment_size, sample_id, queue, semaphore))
        p.start()
        return p

    processes = [start_worker(bed_file) for bed_file in bed_files]

    with open(output_file, 'a') as f:
        while any(p.is_alive() for p in processes) or not queue.empty():
            while not queue.empty():
                bed_file, smoothed_normalized_averages = queue.get()
                line = [sample_id, bed_file.split('.')[0]] + list(map(str, smoothed_normalized_averages))
                f.write("\t".join(line) + "\n")
            for p in processes:
                p.join(0.1)
                if not p.is_alive():
                    processes.remove(p)

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
