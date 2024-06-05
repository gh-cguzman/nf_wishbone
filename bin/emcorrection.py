#!/usr/bin/env python
import argparse
import pysam
import numpy as np
from collections import defaultdict, Counter
import random
import logging
import twobitreader
from scipy.ndimage import gaussian_filter
from multiprocessing import Pool, cpu_count
import os

# Define constants
MIN_FRAGMENT_LENGTH = 50
MAX_FRAGMENT_LENGTH = 550
CHROMOSOMES_TO_PROCESS = set(
    ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
     "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
     "chr18", "chr19", "chr20", "chr21", "chr22"]
)

def parse_args():
    parser = argparse.ArgumentParser(description="Correct end motif bias using 5' end motifs.")
    parser.add_argument("-i", "--input_bam", required=True, help="Input BAM file.")
    parser.add_argument("-r", "--reference_genome", required=True, help="Reference genome file in 2bit format.")
    parser.add_argument("-e", "--exclusion_bed", required=True, help="Exclusion BED file.")
    parser.add_argument("-o", "--output_bam", required=True, help="Output tagged BAM file.")
    parser.add_argument("-n", "--num_simulations", type=int, default=200000000, help="Number of simulations for expected attributes.")
    parser.add_argument("--simulation_rounds", type=int, default=4, help="Number of simulation rounds.")
    parser.add_argument("--output_observed", action='store_true', help="Output TSV file for observed counts.")
    parser.add_argument("--output_expected", action='store_true', help="Output TSV file for expected counts.")
    parser.add_argument("--output_freq", action='store_true', help="Output TSV file for frequency weights.")
    parser.add_argument("--mask_threshold", type=int, default=10, help="Threshold for masking rare attributes.")
    parser.add_argument("--smoothing_sigma", type=float, default=1.0, help="Sigma value for Gaussian smoothing of weights.")
    parser.add_argument("--iqr_factor", type=float, default=1.5, help="Factor for IQR-based outlier limiting.")
    parser.add_argument("--threads", type=int, default=cpu_count(), help="Number of threads to use.")
    parser.add_argument("--sort_bam", action="store_true", help="Sort and merge output BAM file.")
    args = parser.parse_args()
    return args

def setup_logging():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger()
    return logger

def load_exclusion_bed(exclusion_bed_file):
    logger = setup_logging()
    logger.info("Loading exclusion BED file ...")
    excluded_regions = defaultdict(list)
    with open(exclusion_bed_file, "r") as bed:
        for line in bed:
            chrom, start, end = line.strip().split()[:3]
            excluded_regions[chrom].append((int(start), int(end)))
    return excluded_regions

def is_excluded(chrom, start, end, excluded_regions):
    if chrom not in excluded_regions:
        return False
    for region_start, region_end in excluded_regions[chrom]:
        if start < region_end and end > region_start:
            return True
    return False

def count_attributes(bam_file, threads):
    logger = logging.getLogger(__name__)
    logger.info(f"Counting observed attributes from file: {bam_file} ...")
    observed_attributes = Counter()
    with pysam.AlignmentFile(bam_file, "rb", threads=threads) as bam:
        for read in bam:
            if read.is_proper_pair and not read.is_duplicate and MIN_FRAGMENT_LENGTH <= abs(read.template_length) <= MAX_FRAGMENT_LENGTH and read.is_read1:
                seq = read.query_sequence.upper()
                length = abs(read.template_length)
                end_motif = seq[:4]  # 5' end 4 bp motif
                if 'N' not in end_motif:
                    observed_attributes[(length, end_motif)] += 1
    return observed_attributes

def simulate_fragment(ref_genome, fragment_lengths, excluded_regions):
    while True:
        length = random.choice(fragment_lengths)
        chrom = random.choice(list(ref_genome.keys()))
        chrom_len = len(ref_genome[chrom])
        start_pos = random.randint(0, chrom_len - length)
        end_pos = start_pos + length
        if is_excluded(chrom, start_pos, end_pos, excluded_regions):
            continue
        fragment = ref_genome[chrom][start_pos:end_pos].upper()
        end_motif = fragment[:4]
        if 'N' not in end_motif:
            return (length, end_motif)

def worker_simulation(task_params):
    ref_genome, fragment_lengths, excluded_regions, num_simulations, random_seed = task_params
    random.seed(int(random_seed))
    simulated_attributes = Counter()
    for _ in range(num_simulations):
        attribute = simulate_fragment(ref_genome, fragment_lengths, excluded_regions)
        simulated_attributes[attribute] += 1
    return simulated_attributes

def simulate_attributes(ref_genome, fragment_lengths, excluded_regions, num_simulations, threads):
    logger = logging.getLogger(__name__)
    logger.info(f"Simulating {num_simulations} expected attributes ...")

    random_seeds = np.random.randint(0, 1000, size=threads)
    task_params = [(ref_genome, fragment_lengths, excluded_regions, num_simulations // threads, int(seed)) for seed in random_seeds]

    with Pool(threads) as pool:
        results = pool.map(worker_simulation, task_params)

    simulated_attributes = Counter()
    for result in results:
        simulated_attributes.update(result)

    return simulated_attributes

def mask_rare_attributes(observed, threshold=10):
    logger = logging.getLogger(__name__)
    logger.info("Masking rare attributes ...")
    mask = {key: count >= threshold for key, count in observed.items()}
    return mask

def apply_mask(weights, mask):
    return {key: weight if mask.get(key, False) else 1.0 for key, weight in weights.items()}

def smooth_weights(weights, sigma=1.0):
    logger = logging.getLogger(__name__)
    logger.info("Smoothing masked weights ...")
    keys = list(weights.keys())
    lengths = sorted(set(key[0] for key in keys))
    motifs = sorted(set(key[1] for key in keys))
    weight_matrix = np.zeros((len(lengths), len(motifs)))

    length_index = {length: i for i, length in enumerate(lengths)}
    motif_index = {motif: i for i, motif in enumerate(motifs)}

    for (length, motif), weight in weights.items():
        i = length_index[length]
        j = motif_index[motif]
        weight_matrix[i, j] = weight

    smoothed_matrix = gaussian_filter(weight_matrix, sigma=sigma)

    smoothed_weights = {(length, motif): smoothed_matrix[length_index[length], motif_index[motif]] for length, motif in keys}

    return smoothed_weights

def limit_outliers(weights, factor=1.5):
    logger = logging.getLogger(__name__)
    logger.info("Limiting outliers in weights outside 25-75% IQR ...")
    values = np.array(list(weights.values()))
    q1 = np.percentile(values, 25)
    q3 = np.percentile(values, 75)
    iqr = q3 - q1
    upper_limit = q3 + factor * iqr

    return {key: min(value, upper_limit) for key, value in weights.items()}

def compute_weights(observed, expected):
    logger = logging.getLogger(__name__)
    logger.info("Computing bias correction weights ...")
    return {key: expected.get(key, 1.0) / observed[key] if expected.get(key, 0) > 0 else 1.0 for key in observed}

def apply_weights_to_chromosome(input_bam, weights, output_bam_prefix, chromosome, threads):
    logger = logging.getLogger(__name__)
    logger.info(f"Applying weights to chromosome {chromosome} and writing to temporary BAM file ...")
    tmp_output_bam = f"{output_bam_prefix}_{chromosome}.bam"
    with pysam.AlignmentFile(input_bam, "rb", threads=threads) as bam, \
         pysam.AlignmentFile(tmp_output_bam, "wb", header=bam.header) as out_bam:
        for read in bam.fetch(chromosome):
            if read.is_proper_pair and MIN_FRAGMENT_LENGTH <= abs(read.template_length) <= MAX_FRAGMENT_LENGTH:
                seq = read.query_sequence
                length = len(seq)
                end_motif = seq[:4]
                if 'N' not in end_motif:
                    weight = weights.get((length, end_motif), 1.0)
                    read.set_tag("EM", weight)  # Adding weight as a custom tag
                    out_bam.write(read)
    return tmp_output_bam

def load_reference_genome(reference_genome_file):
    logger = logging.getLogger(__name__)
    logger.info("Loading reference genome ...")
    ref_genome = twobitreader.TwoBitFile(reference_genome_file)
    return ref_genome

def sort_bam(temp_bam, threads):
    sorted_temp_file = f"{temp_bam}.sorted.bam"
    pysam.sort("-o", sorted_temp_file, temp_bam, threads=threads)
    os.remove(temp_bam)  # Clean up the temporary unsorted BAM file
    return sorted_temp_file

def sort_and_merge_bams(temp_bam_files, output_bam, threads):
    logger = setup_logging()
    logger.info(f"Sorting temporary BAM files ...")

    # Sort each temporary BAM file in parallel
    with Pool(threads) as pool:
        sorted_temp_files = pool.starmap(sort_bam, [(temp_bam, threads) for temp_bam in temp_bam_files])

    # Merge sorted BAM files
    logger.info(f"Merging temporary BAM files into final BAM file ...")
    pysam.merge("-f", output_bam, *sorted_temp_files, threads=threads)

    # Clean up sorted temporary BAM files
    for sorted_temp_file in sorted_temp_files:
        os.remove(sorted_temp_file)

def write_to_tsv(data, output_file):
    with open(output_file, 'w') as f:
        for key, value in data.items():
            length, motif = key
            f.write(f"{length}\t{motif}\t{value}\n")

def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)
    logger.info(f"Starting end motif bias correction process for BAM {args.input_bam} ...")

    # Load reference genome
    ref_genome = load_reference_genome(args.reference_genome)

    # Load exclusion regions
    excluded_regions = load_exclusion_bed(args.exclusion_bed)

    # Count observed attributes
    observed_attributes = count_attributes(args.input_bam, args.threads)

    # Get fragment lengths for simulation
    fragment_lengths = [length for length, motif in observed_attributes.keys()]

    # Simulate expected attributes
    expected_attributes = simulate_attributes(ref_genome, fragment_lengths, excluded_regions, args.num_simulations, args.threads)

    # Compute bias correction weights
    weights = compute_weights(observed_attributes, expected_attributes)

    # Mask rare or absent attributes
    mask = mask_rare_attributes(observed_attributes, threshold=args.mask_threshold)
    masked_weights = apply_mask(weights, mask)

    # Smooth weights
    smoothed_weights = smooth_weights(masked_weights, sigma=args.smoothing_sigma)

    # Limit outliers in weights
    limited_weights = limit_outliers(smoothed_weights, factor=args.iqr_factor)

    # Apply weights to chromosomes
    bam_file = args.input_bam
    chromosomes = [chrom for chrom in pysam.AlignmentFile(bam_file).references if chrom in CHROMOSOMES_TO_PROCESS]

    # Process each chromosome in parallel
    with Pool(args.threads) as pool:
        temp_bam_files = pool.starmap(apply_weights_to_chromosome, [(bam_file, limited_weights, args.output_bam, chrom, args.threads) for chrom in chromosomes])

    # Sort and merge temporary BAM files if --sort is specified
    if args.sort_bam:
        sort_and_merge_bams(temp_bam_files, args.output_bam, args.threads)
    else:
        logger.info("Skipping sorting and merging of temporary BAM files.")

    # Write observed and expected counts to TSV files if specified
    if args.output_observed:
        write_to_tsv(observed_attributes, 'observed.tsv')
    if args.output_expected:
        write_to_tsv(expected_attributes, 'expected.tsv')

    logger.info("Completed ...")

if __name__ == "__main__":
    main()