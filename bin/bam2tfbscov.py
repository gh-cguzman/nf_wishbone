#!/usr/bin/env python3
import argparse
import logging
import pandas as pd
import pysam
import numpy as np
from pyranges import PyRanges
from scipy.signal import savgol_filter

def setup_logging():
    """Setup logging configuration."""
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def read_bed_file(file):
    """Read a single BED file and return a DataFrame."""
    logging.info(f"Reading BED file: {file}")
    bed = pd.read_csv(file, sep='\t', header=None, skiprows=1)
    bed['source_file'] = file
    return bed

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

def build_pyranges(bed):
    """Build PyRanges object from the BED dataframe."""
    logging.info("Building PyRanges object from BED regions")
    if bed.shape[1] > 4:
        bed = bed.loc[:, [0, 1, 2, 'source_file']]  # Keep only the first four columns if there are more
    bed.columns = ['Chromosome', 'Start', 'End', 'source_file']
    pr = PyRanges(bed)
    return pr

def process_bed_region(bamfile, bed_region, min_fragment_size, max_fragment_size):
    """Process a single BED region and calculate the GC content."""
    chrom = bed_region.Chromosome
    start = bed_region.Start
    end = bed_region.End
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

    return array

def process_bam(bam_file, bed_pr, min_fragment_size, max_fragment_size):
    """Process BAM file and calculate the GC content in extended BED regions."""
    logging.info(f"Processing BAM file: {bam_file}")
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    bed_df = bed_pr.df.copy()
    bed_df['array'] = [np.zeros(end - start + 1) for start, end in zip(bed_df['Start'], bed_df['End'])]

    for bed_region in bed_df.itertuples():
        bed_df.at[bed_region.Index, 'array'] = process_bed_region(bamfile, bed_region, min_fragment_size, max_fragment_size)

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
    result = {}
    for source_file, group in bed_df.groupby('source_file'):
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
    
    combined_averages = {}

    for i in range(0, len(bed_files), 5):
        batch_files = bed_files[i:i+5]
        batch_beds = []
        for bed_file in batch_files:
            bed = read_bed_file(bed_file)
            extended_bed = extend_bed_regions(bed, window=5000)
            filtered_bed = filter_chr(extended_bed)
            sorted_bed = human_sort_bed(filtered_bed)
            batch_beds.append(sorted_bed)
        
        combined_batch_bed = pd.concat(batch_beds, ignore_index=True)
        batch_bed_pr = build_pyranges(combined_batch_bed)
        processed_bed_df = process_bam(bam_file, batch_bed_pr, min_fragment_size, max_fragment_size)
        filtered_bed_df = filter_high_coverage_bins(processed_bed_df, threshold=10)
        batch_averages = average_coverage_by_source(filtered_bed_df)
        combined_averages.update(batch_averages)
    
    smoothed_normalized_averages = smooth_and_normalize(combined_averages)
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
