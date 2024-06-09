#!/usr/bin/env python3
import os
import pysam
import pandas as pd
import numpy as np
from scipy.signal import savgol_filter
from scipy.sparse import lil_matrix
import argparse
import logging
import time
import gc

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Generator to read BED files line by line
def read_bed_file(bed_file):
    logging.info(f'Reading BED file: {bed_file}')
    with open(bed_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                yield {'chr': parts[0], 'start': int(parts[1]), 'end': int(parts[2]), 'bed_file': os.path.basename(bed_file)}

# Function to filter regions to include only chr1-chr22
def filter_chromosomes(data):
    valid_chromosomes = {f'chr{i}' for i in range(1, 23)}
    return (row for row in data if row['chr'] in valid_chromosomes)

# Function to sort regions by human chromosome order
def human_sort(data):
    data = sorted(data, key=lambda row: (int(row['chr'][3:]), row['start'], row['end']))
    for row in data:
        yield row

# Function to extend regions by +/- extension_size bp
def extend_regions(data, extension_size):
    for row in data:
        row['start'] -= extension_size
        row['end'] += extension_size
        yield row

# Modified concatenate_bed_files function to use piped transformations
def concatenate_and_process_bed_files(bed_files, extension_size):
    start_time = time.time()
    data = []
    for bed_file in bed_files:
        try:
            data.extend(read_bed_file(bed_file))
        except Exception as e:
            logging.error(f'Error reading {bed_file}: {e}')
            continue
    
    data = filter_chromosomes(data)
    data = human_sort(data)
    data = list(extend_regions(data, extension_size))
    
    concatenated_df = pd.DataFrame(data, columns=['chr', 'start', 'end', 'bed_file'])
    end_time = time.time()
    logging.info(f'Concatenated and processed BED files in {end_time - start_time:.2f} seconds')
    return concatenated_df

# Generator to process reads from BAM file
def process_bam_file(bam_file, regions):
    logging.info('Processing BAM file')
    start_time = time.time()
    bam = pysam.AlignmentFile(bam_file, 'rb')
    
    coverage_dict = {index: lil_matrix((1, row['end'] - row['start'] + 1), dtype='float32') for index, row in regions.iterrows()}
    
    region_index = 0
    regions_len = len(regions)
    
    regions_start = regions['start'].values
    regions_end = regions['end'].values
    regions_chr = regions['chr'].values
    
    for read in bam.fetch():
        if (read.is_proper_pair and not read.is_duplicate and not read.is_supplementary and 
                110 <= read.template_length <= 230):
            read_chr = read.reference_name
            midpoint = read.reference_start + read.template_length // 2
            
            initial_region_index = region_index
            
            while region_index < regions_len:
                region_chr = regions_chr[region_index]
                region_start = regions_start[region_index]
                region_end = regions_end[region_index]
                
                if read_chr != region_chr:
                    break
                
                if region_start <= midpoint <= region_end:
                    pos_in_array = midpoint - region_start
                    if 0 <= pos_in_array < coverage_dict[region_index].shape[1]:
                        coverage_dict[region_index][0, pos_in_array] += read.get_tag('GC')
                    # Continue to the next region without breaking
                    region_index += 1
                elif midpoint < region_start:
                    break
                elif midpoint > region_end:
                    initial_region_index += 1
                    region_index = initial_region_index
            
            # Reset the region_index to the initial index for the next read
            region_index = initial_region_index

    bam.close()
    end_time = time.time()
    logging.info(f'Processed BAM file in {end_time - start_time:.2f} seconds')
    return coverage_dict

# Function to filter high coverage positions
def filter_high_coverage(coverage_dict):
    start_time = time.time()
    for region_id, coverage_matrix in coverage_dict.items():
        coverage_array = coverage_matrix.toarray().flatten()
        mean = np.mean(coverage_array)
        std = np.std(coverage_array)
        high_coverage_indices = coverage_array > mean + 10 * std
        coverage_matrix[0, high_coverage_indices] = 0
    end_time = time.time()
    logging.info(f'Filtered high coverage in {end_time - start_time:.2f} seconds')
    return coverage_dict

# Function to compute mean coverage profile for each BED file
def compute_mean_profiles(regions, coverage_dict):
    start_time = time.time()
    bed_file_groups = regions.groupby('bed_file')
    mean_profiles = {}
    for bed_file, group in bed_file_groups:
        arrays = [coverage_dict[i].toarray().flatten() for i in group.index]
        stacked_arrays = np.vstack(arrays)  # Stack arrays vertically
        mean_profiles[bed_file] = np.mean(stacked_arrays, axis=0)  # Compute the mean along columns
    end_time = time.time()
    logging.info(f'Computed mean profiles in {end_time - start_time:.2f} seconds')
    return mean_profiles

# Function to combine mean profiles from multiple batches
def combine_mean_profiles(all_mean_profiles):
    combined_mean_profiles = {}

    for mean_profiles in all_mean_profiles:
        for bed_file, profile in mean_profiles.items():
            if bed_file in combined_mean_profiles:
                combined_mean_profiles[bed_file].append(profile)
            else:
                combined_mean_profiles[bed_file] = [profile]
    
    # Average the combined profiles
    for bed_file, profiles in combined_mean_profiles.items():
        stacked_profiles = np.vstack(profiles)
        combined_mean_profiles[bed_file] = np.mean(stacked_profiles, axis=0)
    
    return combined_mean_profiles

# Function to normalize and extract the central region
def normalize_and_extract_center(mean_profiles, extension_size):
    start_time = time.time()
    central_profiles = {}
    center_start = extension_size - 1000
    center_end = extension_size + 1001

    for bed_file, profile in mean_profiles.items():
        mean_profile = np.mean(profile)
        if mean_profile != 0:
            normalized = profile / mean_profile
            central_profiles[bed_file] = normalized[center_start:center_end]
        else:
            central_profiles[bed_file] = np.zeros(center_end - center_start)  # Avoid division by zero
    end_time = time.time()
    logging.info(f'Normalized and extracted central regions in {end_time - start_time:.2f} seconds')
    return central_profiles

# Function to smooth coverage profiles
def smooth_profiles(central_profiles):
    start_time = time.time()
    smoothed_profiles = {}
    for bed_file, profile in central_profiles.items():
        smoothed = savgol_filter(profile, window_length=165, polyorder=3)
        smoothed_profiles[bed_file] = smoothed
    end_time = time.time()
    logging.info(f'Smoothed profiles in {end_time - start_time:.2f} seconds')
    return smoothed_profiles

# Function to output the +/- 1000bp region to a TSV file
def output_profiles(smoothed_profiles, output_file, sample_id):
    start_time = time.time()
    with open(output_file, 'w') as f:
        f.write('Sample_ID\tBED_File\t' + '\t'.join(map(str, range(-1000, 1001))) + '\n')
        for bed_file, profile in smoothed_profiles.items():
            f.write(f'{sample_id}\t{bed_file.split(".")[0]}\t' + '\t'.join(map(str, profile)) + '\n')
    end_time = time.time()
    logging.info(f'Output profiles to TSV file in {end_time - start_time:.2f} seconds')

def batch_bed_files(bed_files, batch_size=30):
    """Yield successive n-sized chunks from list of bed_files."""
    for i in range(0, len(bed_files), batch_size):
        yield bed_files[i:i + batch_size]

# Main function to execute the workflow
def main(bam_file, bed_files, sample_id, output_file, extension_size=5000, batch_size=30):
    start_time = time.time()
    
    all_mean_profiles = []

    for batch_num, bed_file_batch in enumerate(batch_bed_files(bed_files, batch_size)):
        logging.info(f'Processing batch {batch_num + 1}')
        
        extended_bed = concatenate_and_process_bed_files(bed_file_batch, extension_size)
        
        coverage_dict = process_bam_file(bam_file, extended_bed)
        #gc.collect()  # Collect garbage to free up memory
        
        filtered_coverage_dict = filter_high_coverage(coverage_dict)
        #gc.collect()  # Collect garbage to free up memory
        
        mean_profiles = compute_mean_profiles(extended_bed, filtered_coverage_dict)
        all_mean_profiles.append(mean_profiles)

        # Delete intermediate variables to free memory
        del extended_bed
        del coverage_dict
        del filtered_coverage_dict
        gc.collect()

    # Combine mean profiles from all batches
    combined_mean_profiles = combine_mean_profiles(all_mean_profiles)
    central_profiles = normalize_and_extract_center(combined_mean_profiles, extension_size)
    smoothed_profiles = smooth_profiles(central_profiles)
    
    output_profiles(smoothed_profiles, output_file, sample_id)
    
    end_time = time.time()
    logging.info(f'Total execution time: {end_time - start_time:.2f} seconds')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process BAM and BED files to compute coverage profiles.')
    parser.add_argument('-B', '--bam-file', required=True, help='Input BAM file.')
    parser.add_argument('-b', '--bed-files', nargs='+', required=True, help='Input BED files.')
    parser.add_argument('-s', '--sample-id', required=True, help='Sample ID for the output file.')
    parser.add_argument('-o', '--cov-output', required=True, help='Output file for the coverage profiles.')
    parser.add_argument('-e', '--extension-size', type=int, default=5000, help='Extension size in bp (default: 5000)')
    parser.add_argument('-bs', '--batch-size', type=int, default=30, help='Batch size for processing BED files (default: 30)')
    
    args = parser.parse_args()
    main(args.bam_file, args.bed_files, args.sample_id, args.cov_output, args.extension_size, args.batch_size)
