#!/usr/bin/env python3
import argparse
import logging
import pandas as pd
import pysam
import numpy as np
from intervaltree import Interval, IntervalTree
from scipy.signal import savgol_filter

def setup_logging():
    """Setup logging configuration."""
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def concat_bed_files(bed_files):
    """Concatenate multiple BED files and keep track of the source file."""
    bed_frames = []
    for file in bed_files:
        logging.info(f"Reading BED file: {file}")
        bed = pd.read_csv(file, sep='\t', header=None, skiprows=1)
        bed['source_file'] = file
        bed_frames.append(bed)
    concatenated = pd.concat(bed_frames, ignore_index=True)
    return concatenated

def extend_bed_regions(bed, window=5000):
    """Extend BED regions by ±window bp."""
    logging.info(f"Extending BED regions by ±{window} bp")
    bed[1] = bed[1] - window
    bed[2] = bed[2] + window
    bed[1] = bed[1].clip(lower=0)  # Ensure start is not less than 0
    return bed

def filter_chr(bed):
    """Keep only regions in chr1 to chr22."""
    logging.info("Filtering BED regions to keep only chr1 to chr22")
    chr_filter = ['chr' + str(i) for i in range(1, 23)]
    bed_filtered = bed[bed[0].isin(chr_filter)].copy()
    return bed_filtered

def human_sort_bed(bed):
    """Sort BED file in human-readable order (chr1, chr2, ..., chr10, chr11, ..., chr22)."""
    logging.info("Sorting BED regions in human-readable order")
    bed.loc[:, 'chrom'] = pd.Categorical(bed[0], categories=[f'chr{i}' for i in range(1, 23)], ordered=True)
    sorted_bed = bed.sort_values(['chrom', 1, 2])
    return sorted_bed.drop(columns='chrom')

def build_interval_tree(bed):
    """Build an interval tree from the BED dataframe."""
    logging.info("Building interval tree from BED regions")
    interval_trees = {f'chr{i}': IntervalTree() for i in range(1, 23)}
    for index, row in bed.iterrows():
        chrom = row[0]
        start = row[1]
        end = row[2]
        interval_trees[chrom][start:end] = index
    return interval_trees

def process_bam(bam_file, bed, interval_trees, min_fragment_size=100, max_fragment_size=200):
    """Process BAM file and calculate the GC content in extended BED regions."""
    logging.info(f"Processing BAM file: {bam_file}")
    bed['array'] = bed.apply(lambda x: np.zeros(x[2] - x[1] + 1), axis=1)

    bamfile = pysam.AlignmentFile(bam_file, "rb")
    read_count = 0
    
    for read in bamfile.fetch():
        read_count += 1
        if read_count % 10_000_000 == 0:
            logging.info(f"Processed {read_count // 1_000_000}M reads")
        
        if not read.is_proper_pair or read.is_duplicate or read.is_supplementary:
            continue
        
        if min_fragment_size <= read.template_length <= max_fragment_size:
            midpoint = (read.reference_start + read.template_length) // 2
            gc_value = 1 / read.get_tag('GC')

            if read.reference_name in interval_trees:
                overlaps = interval_trees[read.reference_name][midpoint]
                for interval in overlaps:
                    index = interval.data
                    start = bed.at[index, 1]
                    pos_in_array = midpoint - start
                    bed.at[index, 'array'][pos_in_array] += gc_value
    
    bamfile.close()
    logging.info(f"Finished processing BAM file with {read_count} reads")
    return bed

def filter_high_coverage_bins(bed, threshold=10):
    """Filter out regions with extremely high coverage (above threshold standard deviations)."""
    logging.info(f"Filtering out regions with coverage above {threshold} standard deviations")
    all_values = np.concatenate(bed['array'].values)
    mean = np.mean(all_values)
    std = np.std(all_values)
    cutoff = mean + threshold * std
    
    bed['filtered_array'] = np.empty((len(bed),), dtype=object)
    for index, row in bed.iterrows():
        array = row['array']
        array[array > cutoff] = 0
        bed.at[index, 'filtered_array'] = array
    return bed

def average_coverage_by_source(bed):
    """Compute the average coverage for each source file."""
    logging.info("Computing average coverage for each source file")
    result = {}
    for source_file, group in bed.groupby('source_file'):
        filtered_arrays = np.vstack(group['filtered_array'].values)
        avg_coverage = np.mean(filtered_arrays, axis=0)
        result[source_file] = avg_coverage
    return result

def smooth_and_normalize(averages):
    """Smooth and normalize arrays."""
    logging.info("Smoothing and normalizing coverage")
    smoothed_averages = {}
    for source_file, avg_coverage in averages.items():
        # Normalize to a mean coverage of 1 across the -5000 to +5000 region
        mean_coverage = np.mean(avg_coverage)
        if mean_coverage > 0:
            normalized_coverage = avg_coverage / mean_coverage
        else:
            normalized_coverage = np.zeros_like(avg_coverage)
        # Extract -1000 to +1000 region
        coverage_region = normalized_coverage[4000:6001]
        # Smooth the coverage region
        filtered_array = savgol_filter(coverage_region, window_length=165, polyorder=3)
        smoothed_averages[source_file] = filtered_array
    return smoothed_averages

def save_averaged_coverage_to_tsv(averages, output_file, sample_id):
    """Save the averaged coverage data to a TSV file."""
    logging.info(f"Saving averaged coverage to TSV file: {output_file}")
    with open(output_file, 'w') as f:
        header = ["Sample_ID", "BED_File"] + [str(i) for i in range(-1000, 1001)]
        f.write("\t".join(header) + "\n")
        
        for source_file, avg_coverage in averages.items():
            bed_file = source_file.split('.')[0]
            line = [sample_id, bed_file] + list(map(str, avg_coverage))
            f.write("\t".join(line) + "\n")

def main(args):
    setup_logging()
    
    bed_files = args.bed_files
    bam_file = args.bam_file
    output_file = args.output_file
    sample_id = args.sample_id
    min_fragment_size = args.min_fragment_size
    max_fragment_size = args.max_fragment_size
    
    concatenated_bed = concat_bed_files(bed_files)
    extended_bed = extend_bed_regions(concatenated_bed, window=5000)
    filtered_bed = filter_chr(extended_bed)
    sorted_bed = human_sort_bed(filtered_bed)
    
    interval_trees = build_interval_tree(sorted_bed)
    processed_bed = process_bam(bam_file, sorted_bed, interval_trees, min_fragment_size, max_fragment_size)
    
    filtered_bed = filter_high_coverage_bins(processed_bed, threshold=10)
    
    averaged_coverage = average_coverage_by_source(filtered_bed)
    
    smoothed_normalized_averages = smooth_and_normalize(averaged_coverage)
    
    save_averaged_coverage_to_tsv(smoothed_normalized_averages, output_file, sample_id)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process BAM and BED files.")
    parser.add_argument("-b", "--bed-files", nargs='+', required=True, help="List of BED files.")
    parser.add_argument("-B", "--bam-file", required=True, help="Input BAM file.")
    parser.add_argument("-o", "--output-file", required=True, help="Output TSV file.")
    parser.add_argument("-s", "--sample-id", required=True, help="Sample ID.")
    parser.add_argument("-m", "--min-fragment-size", type=int, default=110, help="Minimum fragment size.")
    parser.add_argument("-M", "--max-fragment-size", type=int, default=230, help="Maximum fragment size.")
    
    args = parser.parse_args()
    main(args)
