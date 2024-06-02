#!/usr/bin/env python
import argparse
import pysam
import numpy as np
from collections import defaultdict
from multiprocessing import Pool
from functools import partial
import random
import logging
import twobitreader
from scipy.ndimage import gaussian_filter

# Define constants
MIN_FRAGMENT_LENGTH = 50
MAX_FRAGMENT_LENGTH = 550
RANDOM_SEED = random.randint(1, 999)

def parse_args():
    parser = argparse.ArgumentParser(description="Correct end motif bias using 5' end motifs.")
    parser.add_argument("-i", "--input_bam", required=True, help="Input BAM file.")
    parser.add_argument("-r", "--reference_genome", required=True, help="Reference genome file in 2bit format.")
    parser.add_argument("-e", "--exclusion_bed", required=True, help="Exclusion BED file.")
    parser.add_argument("-o", "--output_bam", help="Output tagged BAM file.")
    parser.add_argument("-n", "--num_simulations", type=int, default=200000000, help="Number of simulations for expected attributes.")
    parser.add_argument("--simulation_rounds", type=int, default=4, help="Number of simulation rounds.")
    parser.add_argument("--output_observed", action='store_true', help="Output TSV file for observed counts.")
    parser.add_argument("--output_expected", action='store_true', help="Output TSV file for expected counts.")
    parser.add_argument("--output_freq", action='store_true', help="Output TSV file for frequency weights.")
    parser.add_argument("--mask_threshold", type=int, default=10, help="Threshold for masking rare attributes.")
    parser.add_argument("--smoothing_sigma", type=float, default=1.0, help="Sigma value for Gaussian smoothing of weights.")
    parser.add_argument("--iqr_factor", type=float, default=1.5, help="Factor for IQR-based outlier limiting.")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use.")
    args = parser.parse_args()
    return args

def setup_logging():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger()
    return logger

def count_attributes(bam_file, max_reads=50000000):
    logger = setup_logging()
    observed_attributes = defaultdict(int)
    logger.info(f"Counting fragment attributes from BAM file, limited to {max_reads} reads ...")
    read_count = 0

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            if read.is_proper_pair and not read.is_duplicate and MIN_FRAGMENT_LENGTH <= abs(read.template_length) <= MAX_FRAGMENT_LENGTH and read.is_read1:
                seq = read.query_sequence.upper()
                length = abs(read.template_length)
                end_motif = seq[:4]  # 5' end 4 bp motif
                if 'N' not in end_motif:
                    observed_attributes[(length, end_motif)] += 1
                    read_count += 1
                    if read_count >= max_reads:
                        logger.info(f"Reached the limit of {max_reads} reads. Continuing ...")
                        break

    return observed_attributes

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

def simulate_fragment(ref_genome, fragment_lengths, excluded_regions, _):
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

def simulate_attributes(ref_genome, fragment_lengths, excluded_regions, num_simulations, rounds, threads):
    logger = setup_logging()
    simulated_attributes = defaultdict(int)
    logger.info("Simulating expected attributes ...")
    for round in range(rounds):
        logger.info(f"Simulation round {round + 1} of {rounds}: Simulating {len(range(num_simulations // rounds))} fragments ...")
        pool = Pool(threads)
        simulate_partial = partial(simulate_fragment, ref_genome, fragment_lengths, excluded_regions)
        results = pool.map(simulate_partial, range(num_simulations // rounds))
        pool.close()
        pool.join()
        for length, end_motif in results:
            simulated_attributes[(length, end_motif)] += 1
    return simulated_attributes

def mask_rare_attributes(observed, threshold=10):
    logger = setup_logging()
    logger.info("Masking rare attributes ...")
    mask = {}
    for key, count in observed.items():
        mask[key] = count >= threshold
    return mask

def apply_mask(weights, mask):
    logger = setup_logging()
    logger.info("Applying binary mask ...")
    masked_weights = {}
    for key, weight in weights.items():
        if mask.get(key, False):
            masked_weights[key] = weight
        else:
            masked_weights[key] = 1.0  # Default weight for masked attributes
    return masked_weights

def smooth_weights(weights, sigma=1.0):
    logger = setup_logging()
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

    smoothed_weights = {}
    for (length, motif), weight in weights.items():
        i = length_index[length]
        j = motif_index[motif]
        smoothed_weights[(length, motif)] = smoothed_matrix[i, j]

    return smoothed_weights

def limit_outliers(weights, factor=1.5):
    logger = setup_logging()
    logger.info("Limiting outliers in weights outside 25-75% IQR ...")
    values = np.array(list(weights.values()))
    q1 = np.percentile(values, 25)
    q3 = np.percentile(values, 75)
    iqr = q3 - q1
    upper_limit = q3 + factor * iqr

    limited_weights = {k: min(v, upper_limit) for k, v in weights.items()}
    return limited_weights

def compute_weights(observed, expected):
    logger = setup_logging()
    weights = {}
    logger.info("Computing bias correction weights ...")
    for key in observed:
        if key in expected and expected[key] > 0:
            weights[key] = expected[key] / observed[key]
        else:
            weights[key] = 1.0  # Default weight for rare/absent attributes
    return weights

def apply_weights(bam_file, weights, output_bam):
    logger = setup_logging()
    logger.info("Applying weights and writing to output BAM file ...")
    with pysam.AlignmentFile(bam_file, "rb") as bam, \
         pysam.AlignmentFile(output_bam, "wb", header=bam.header) as out_bam:
        for read in bam:
            if read.is_proper_pair and MIN_FRAGMENT_LENGTH <= abs(read.template_length) <= MAX_FRAGMENT_LENGTH:
                seq = read.query_sequence
                length = len(seq)
                end_motif = seq[:4]
                if 'N' not in end_motif:
                    weight = weights.get((length, end_motif), 1.0)
                    read.set_tag("EM", weight)  # Adding weight as a custom tag
                    out_bam.write(read)

def load_reference_genome(reference_genome_file):
    logger = setup_logging()
    logger.info("Loading reference genome ...")
    ref_genome = twobitreader.TwoBitFile(reference_genome_file)
    return ref_genome

def write_to_tsv(data, output_file):
    with open(output_file, 'w') as f:
        for key, value in data.items():
            length, motif = key
            f.write(f"{length}\t{motif}\t{value}\n")

def main():
    args = parse_args()
    logger = setup_logging()
    logger.info(f"Starting end motif bias correction process for BAM {args.input_bam} ...")

    # Load reference genome
    ref_genome = load_reference_genome(args.reference_genome)

    # Load exclusion regions
    excluded_regions = load_exclusion_bed(args.exclusion_bed)

    # Count observed attributes
    observed_attributes = count_attributes(args.input_bam)

    # Get fragment lengths for simulation
    fragment_lengths = [length for (length, motif) in observed_attributes.keys()]

    # Simulate expected attributes
    expected_attributes = simulate_attributes(ref_genome, fragment_lengths, excluded_regions, args.num_simulations, args.simulation_rounds, args.threads)

    # Compute bias correction weights
    weights = compute_weights(observed_attributes, expected_attributes)

    # Mask rare or absent attributes
    mask = mask_rare_attributes(observed_attributes, threshold=args.mask_threshold)
    masked_weights = apply_mask(weights, mask)

    # Smooth weights
    smoothed_weights = smooth_weights(masked_weights, sigma=args.smoothing_sigma)

    # Limit outliers in weights
    limited_weights = limit_outliers(smoothed_weights, factor=args.iqr_factor)

    # Apply weights and write to output BAM
    if args.output_bam:
        apply_weights(args.input_bam, limited_weights, args.output_bam)

    # Write smoothed weights to TSV by default
    if args.output_freq:
        write_to_tsv(limited_weights, 'smoothed_weights.tsv')

    # Write observed and expected counts to TSV files if specified
    if args.output_observed:
        write_to_tsv(observed_attributes, 'observed.tsv')
    if args.output_expected:
        write_to_tsv(expected_attributes, 'expected.tsv')

    logger.info("Completed ...")

if __name__ == "__main__":
    main()