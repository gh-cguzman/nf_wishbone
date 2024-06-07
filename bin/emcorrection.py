#!/usr/bin/env python3
import collections
import pysam
import multiprocessing
from functools import partial
import argparse
from collections import defaultdict
import logging
import py2bit
import random
from intervaltree import Interval, IntervalTree
import pandas as pd
import numpy as np
from statsmodels.api import GLM, families
from scipy.ndimage import gaussian_filter
import dill as pickle
import os
import shutil

rc = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

def defaultdict_int():
    return defaultdict(int)

def process_read(read, fragment_min, fragment_max):
    if (read.is_read1 and
        read.is_proper_pair and
        not read.is_duplicate and
        fragment_min <= read.template_length <= fragment_max and
        read.reference_name in [f'chr{i}' for i in range(1, 23)]):
        
        fragment_length = read.template_length
        seq = read.query_sequence
        end_motif = None  # Initialize the end_motif variable
        
        if read.flag in (99, 163):
            end_motif = seq[:4]
        elif read.flag in (83, 147):
            end_motif = ''.join(rc[base] for base in reversed(seq[-4:]))

        if end_motif and 'N' not in end_motif:
            return (fragment_length, end_motif)
    
    return None

def process_chromosome(chromosome, bam_file, fragment_min, fragment_max, queue):
    bam = pysam.AlignmentFile(bam_file, "rb")
    results = defaultdict(defaultdict_int)
    fragment_length_counts = defaultdict_int()
    
    for read in bam.fetch(chromosome):
        result = process_read(read, fragment_min, fragment_max)
        if result:
            fragment_length, end_motif = result
            results[fragment_length][end_motif] += 1
            fragment_length_counts[fragment_length] += 1
    
    queue.put((results, fragment_length_counts))

def merge_results(results_list):
    merged = defaultdict(defaultdict_int)
    for results in results_list:
        for fragment_length, motifs in results.items():
            for motif, count in motifs.items():
                merged[fragment_length][motif] += count
    return merged

def merge_counts(counts_list):
    merged = defaultdict_int()
    for counts in counts_list:
        for fragment_length, number in counts.items():
            merged[fragment_length] += number
    return merged

def parallel_process_bam(bam_file, chromosomes, fragment_min, fragment_max, num_cores):
    processes = []
    queue = multiprocessing.Queue()

    for chrom in chromosomes:
        process = multiprocessing.Process(target=process_chromosome, args=(chrom, bam_file, fragment_min, fragment_max, queue))
        processes.append(process)
    
    for process in processes:
        process.start()
    
    results_list = []
    fragment_length_counts_list = []
    for _ in chromosomes:
        results, fragment_length_counts = queue.get()
        results_list.append(results)
        fragment_length_counts_list.append(fragment_length_counts)
    
    for process in processes:
        process.join()

    motif_results = merge_results(results_list)
    fragment_length_counts = merge_counts(fragment_length_counts_list)
    
    return motif_results, fragment_length_counts

def compute_allowed_intervals(chr_length, blacklist_tree):
    allowed_intervals = []
    start = 0
    for interval in sorted(blacklist_tree):
        if start < interval.begin:
            allowed_intervals.append((start, interval.begin))
        start = interval.end
    if start < chr_length:
        allowed_intervals.append((start, chr_length))
    return allowed_intervals

def sample_from_allowed_intervals(allowed_intervals, fragment_length):
    while True:
        interval = random.choice(allowed_intervals)
        if interval[1] - interval[0] >= fragment_length:  # Ensure the interval is large enough
            start = random.randint(interval[0], interval[1] - fragment_length)
            if start + fragment_length <= interval[1]:
                return start

def simulate_motifs(chromosome, twobit_file, fragment_length_counts, observed_motifs, num_simulations, allowed_intervals, queue):
    twobit = py2bit.open(twobit_file)
    results = defaultdict(defaultdict_int)
    fragment_length_counts_sim = defaultdict_int()
    chr_length = twobit.chroms(chromosome)
    
    fragment_lengths = list(fragment_length_counts.keys())
    fragment_length_weights = list(fragment_length_counts.values())

    for _ in range(num_simulations):
        fragment_length = random.choices(fragment_lengths, weights=fragment_length_weights, k=1)[0]
        start = sample_from_allowed_intervals(allowed_intervals, fragment_length)
        
        sequence = twobit.sequence(chromosome, start, start + 4)
        
        if 'N' not in sequence and sequence in observed_motifs[fragment_length]:
            results[fragment_length][sequence] += 1
            fragment_length_counts_sim[fragment_length] += 1

    twobit.close()
    
    queue.put((results, fragment_length_counts_sim))

def parallel_simulate_motifs(twobit_file, chromosomes, fragment_length_counts, observed_motifs, num_simulations, num_cores, allowed_intervals_dict):
    processes = []
    queue = multiprocessing.Queue()
    
    for chrom in chromosomes:
        process = multiprocessing.Process(target=simulate_motifs, args=(chrom, twobit_file, fragment_length_counts, observed_motifs, num_simulations // len(chromosomes), allowed_intervals_dict[chrom], queue))
        processes.append(process)
    
    for process in processes:
        process.start()
    
    results_list = []
    fragment_length_counts_list = []
    for _ in chromosomes:
        results, fragment_length_counts_sim = queue.get()
        results_list.append(results)
        fragment_length_counts_list.append(fragment_length_counts_sim)
    
    for process in processes:
        process.join()

    motif_results = merge_results(results_list)
    fragment_length_counts = merge_counts(fragment_length_counts_list)

    return motif_results, fragment_length_counts

def load_blacklist(bed_file):
    blacklist_trees = defaultdict(IntervalTree)
    with open(bed_file, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split()
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                blacklist_trees[chrom].add(Interval(start, end))
    return blacklist_trees

def precompute_allowed_intervals(chromosomes, twobit_file, blacklist_trees):
    allowed_intervals_dict = {}
    twobit = py2bit.open(twobit_file)
    for chrom in chromosomes:
        chr_length = twobit.chroms(chrom)
        blacklist_tree = blacklist_trees.get(chrom, IntervalTree())
        allowed_intervals_dict[chrom] = compute_allowed_intervals(chr_length, blacklist_tree)
    twobit.close()
    return allowed_intervals_dict

def calculate_simple_weights(observed, expected):
    weights = collections.defaultdict(partial(collections.defaultdict, float))
    for frag_length in observed:
        for motif in observed[frag_length]:
            obs_count = observed[frag_length][motif]
            exp_count = expected[frag_length].get(motif, 0)
            if obs_count <= 20 or exp_count <= 20:
                weights[frag_length][motif] = 1
            else:
                weights[frag_length][motif] = 1 / (exp_count / obs_count)
    return weights

def calculate_simple_fragment_weights(observed, expected):
    weights = defaultdict(float)
    for frag_length in observed:
        obs_count = observed[frag_length]
        exp_count = expected.get(frag_length, 0)
        if obs_count <= 20 or exp_count <= 20:
            weights[frag_length] = 1
        else:
            weights[frag_length] = 1 / (exp_count / obs_count)
    return weights

def calculate_glm_weights(observed, expected):
    weights = collections.defaultdict(partial(collections.defaultdict, float))
    for frag_length in observed:
        data = []
        for motif in observed[frag_length]:
            obs_count = observed[frag_length][motif]
            exp_count = expected[frag_length].get(motif, 0)
            if obs_count > 0 and exp_count > 0:
                data.append((obs_count, exp_count))
        if data:
            obs_counts, exp_counts = zip(*data)
            glm = GLM(obs_counts, exp_counts, family=families.Poisson()).fit()
            for motif in observed[frag_length]:
                weights[frag_length][motif] = 1 / (glm.predict(expected[frag_length].get(motif, 0))[0])
    return weights

def calculate_glm_fragment_weights(observed, expected):
    data = []
    for frag_length in observed:
        obs_count = observed[frag_length]
        exp_count = expected.get(frag_length, 0)
        if obs_count > 0 and exp_count > 0:
            data.append((obs_count, exp_count))
    weights = defaultdict(float)
    if data:
        obs_counts, exp_counts = zip(*data)
        glm = GLM(obs_counts, exp_counts, family=families.Poisson()).fit()
        for frag_length in observed:
            weights[frag_length] = 1 / (glm.predict(expected.get(frag_length, 0))[0])
    return weights

def write_weights_to_tsv(weights, filename):
    with open(filename, 'w') as f:
        f.write("FragmentLength\tMotif\tWeight\n")
        for frag_length in sorted(weights.keys()):
            for motif, weight in sorted(weights[frag_length].items()):
                f.write(f"{frag_length}\t{motif}\t{weight}\n")

def write_fragment_weights_to_tsv(weights, filename):
    with open(filename, 'w') as f:
        f.write("FragmentLength\tWeight\n")
        for frag_length, weight in sorted(weights.items()):
            f.write(f"{frag_length}\t{weight}\n")

def limit_extreme_weights(weights, sod):
    fs = 10 - sod
    non_default_weights = [weight for frag_length in weights for weight in weights[frag_length].values() if weight != 1]
    
    if not non_default_weights:
        return weights
    
    WQ3 = np.percentile(non_default_weights, 75)
    IQRW = np.percentile(non_default_weights, 75) - np.percentile(non_default_weights, 25)
    Wthreshold = WQ3 + fs * IQRW

    for frag_length in weights:
        for motif in weights[frag_length]:
            if weights[frag_length][motif] > Wthreshold:
                weights[frag_length][motif] = Wthreshold
    
    return weights

def limit_extreme_fragment_weights(weights, sod):
    fs = 10 - sod
    non_default_weights = [weight for weight in weights.values() if weight != 1]
    
    if not non_default_weights:
        return weights
    
    WQ3 = np.percentile(non_default_weights, 75)
    IQRW = np.percentile(non_default_weights, 75) - np.percentile(non_default_weights, 25)
    Wthreshold = WQ3 + fs * IQRW

    for frag_length in weights:
        if weights[frag_length] > Wthreshold:
            weights[frag_length] = Wthreshold
    
    return weights

def smooth_weights(weights, sigma):

    motif_to_index, index_to_motif = create_motif_index_mapping(weights)
    max_frag_length = max(weights.keys())
    max_motif_index = max(motif_to_index.values())
    weight_matrix = np.zeros((max_frag_length + 1, max_motif_index + 1))

    for frag_length in weights:
        for motif in weights[frag_length]:
            weight_matrix[frag_length, motif_to_index[motif]] = weights[frag_length][motif]
    
    smoothed_matrix = gaussian_filter(weight_matrix, sigma=sigma)

    smoothed_weights = collections.defaultdict(partial(collections.defaultdict, float))
    for frag_length in weights:
        for motif in weights[frag_length]:
            smoothed_weights[frag_length][motif] = smoothed_matrix[frag_length, motif_to_index[motif]]
    
    return smoothed_weights

def smooth_fragment_weights(weights, sigma):

    max_frag_length = max(weights.keys())
    weight_array = np.zeros(max_frag_length + 1)

    for frag_length in weights:
        weight_array[frag_length] = weights[frag_length]
    
    smoothed_array = gaussian_filter(weight_array, sigma=sigma)

    smoothed_weights = defaultdict(float)
    for frag_length in weights:
        smoothed_weights[frag_length] = smoothed_array[frag_length]
    
    return smoothed_weights

def create_motif_index_mapping(weights):
    motifs = set()
    for frag_length in weights:
        motifs.update(weights[frag_length].keys())
    motif_to_index = {motif: i for i, motif in enumerate(sorted(motifs))}
    index_to_motif = {i: motif for motif, i in motif_to_index.items()}
    return motif_to_index, index_to_motif

def add_tags_to_read(read, motif_weights, fragment_weights):
    fragment_length = read.template_length
    if read.flag in (99, 163):
        end_motif = read.query_sequence[:4]
    elif read.flag in (83, 147):
        end_motif = ''.join(rc[base] for base in reversed(read.query_sequence[-4:]))
    else:
        end_motif = "NNNN"  # Handle unexpected cases
    
    em_weight = motif_weights.get(fragment_length, {}).get(end_motif, 1.0)
    fl_weight = fragment_weights.get(fragment_length, 1.0)
    read.set_tag("EM", em_weight)
    read.set_tag("FL", fl_weight)
    return read

def process_chromosome_with_tags(chromosome, bam_file, motif_weights, fragment_weights, fragment_min, fragment_max, tmp_dir, conn):
    logging.info(f"Processing chromosome: {chromosome}")
    bam_in = pysam.AlignmentFile(bam_file, "rb")
    tagged_bam_file = os.path.join(tmp_dir, f"{chromosome}_tagged.bam")
    bam_out = pysam.AlignmentFile(tagged_bam_file, "wb", template=bam_in)

    for read in bam_in.fetch(chromosome):
        if (read.is_proper_pair and
            not read.is_supplementary and
            read.reference_name in [f'chr{i}' for i in range(1, 23)] and
            fragment_min <= read.template_length <= fragment_max and
            'N' not in read.query_sequence[:4]):
            
            read = add_tags_to_read(read, motif_weights, fragment_weights)
            bam_out.write(read)

    bam_in.close()
    bam_out.close()
    conn.send(tagged_bam_file)
    conn.send(None)  # Signal that this chromosome is done
    conn.close()

def sort_and_merge_bam_files(output_bam_file, tmp_dir, num_chromosomes, pipes, num_cores):
    logging.info(f"Starting to sort and merge BAM files into {output_bam_file}")
    sorted_bam_files = []

    finished_chromosomes = 0
    while finished_chromosomes < num_chromosomes:
        for pipe in pipes:
            tagged_bam_file = pipe.recv()
            if tagged_bam_file is None:
                finished_chromosomes += 1
                logging.info(f"Finished processing chromosome: {finished_chromosomes}/{num_chromosomes}")
                pipes.remove(pipe)
            else:
                # Sort the temporary BAM file
                sorted_bam_file = os.path.join(tmp_dir, f"sorted_{os.path.basename(tagged_bam_file)}")
                pysam.sort("-@", str(num_cores), "-o", sorted_bam_file, tagged_bam_file)
                sorted_bam_files.append(sorted_bam_file)
    
    if sorted_bam_files:
        # Merge the sorted BAM files
        tmp_merged_bam_file = os.path.join(tmp_dir, "merged_sorted.bam")
        pysam.merge("-@", str(num_cores), tmp_merged_bam_file, *sorted_bam_files)
        os.rename(tmp_merged_bam_file, output_bam_file)
        
        for bam_file in sorted_bam_files:
            os.remove(bam_file)  # Remove temporary sorted BAM files
    
    logging.info(f"Finished sorting and merging BAM files into {output_bam_file}")

def parallel_process_bam_with_tags(bam_file, output_bam_file, chromosomes, motif_weights, fragment_weights, fragment_min, fragment_max, num_cores, tmp_dir):
    processes = []
    parent_conns = []
    
    for chromosome in chromosomes:
        parent_conn, child_conn = multiprocessing.Pipe()
        p = multiprocessing.Process(target=process_chromosome_with_tags, args=(chromosome, bam_file, motif_weights, fragment_weights, fragment_min, fragment_max, tmp_dir, child_conn))
        processes.append(p)
        parent_conns.append(parent_conn)
        p.start()

    merger_process = multiprocessing.Process(target=sort_and_merge_bam_files, args=(output_bam_file, tmp_dir, len(chromosomes), parent_conns, num_cores))
    merger_process.start()

    for p in processes:
        p.join()

    merger_process.join()

def main():
    parser = argparse.ArgumentParser(description='Process BAM file for 5\'-end motif bias correction and simulate end motifs from genome.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input BAM file')
    parser.add_argument('-o', '--output-observed', type=str, help='Output file for observed results')
    parser.add_argument('-e', '--output-expected', type=str, help='Output file for expected (simulated) results')
    parser.add_argument('-fo', '--output-fragment-observed', type=str, help='Output file for observed fragment length counts')
    parser.add_argument('-fe', '--output-fragment-expected', type=str, help='Output file for expected (simulated) fragment length counts')
    parser.add_argument('-cw', '--correction-weights', type=str, required=True, help='Output file for correction weights')
    parser.add_argument('-fcw', '--fragment-correction-weights', type=str, required=True, help='Output file for fragment length correction weights')
    parser.add_argument('-cm', '--correction-method', type=str, choices=['simple', 'glm'], default='simple', help='Method for calculating correction weights (default: simple)')
    parser.add_argument('-f', '--fragment-min', type=int, default=50, help='Minimum fragment length (default: 50)')
    parser.add_argument('-F', '--fragment-max', type=int, default=550, help='Maximum fragment length (default: 550)')
    parser.add_argument('-c', '--chromosomes', type=str, nargs='+', default=[f'chr{i}' for i in range(1, 23)], help='List of chromosomes to process (default: chr1 to chr22)')
    parser.add_argument('-l', '--log', type=str, default="process.log", help='Log file (default: process.log)')
    parser.add_argument('-n', '--num-cores', type=int, default=multiprocessing.cpu_count(), help='Number of CPU cores to use (default: all available cores)')
    parser.add_argument('-t', '--twobit', type=str, required=True, help='2bit file of the reference genome')
    parser.add_argument('-s', '--num-simulations', type=int, default=200000000, help='Number of simulations to perform (default: 200M)')
    parser.add_argument('-b', '--blacklist', type=str, required=True, help='BED file with blacklisted regions')
    parser.add_argument('-sod', '--stringency', type=int, default=3, help='Outlier detection stringency (default: 3)')
    parser.add_argument('-sg', '--smoothing-sigma', type=float, default=2.0, help='Sigma for Gaussian smoothing (default: 2.0)')
    parser.add_argument('-obam', '--output-bam', type=str, help='Output BAM file with added tags')

    args = parser.parse_args()

    logging.basicConfig(filename=args.log, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    bam_file = args.input
    output_observed = args.output_observed
    output_expected = args.output_expected
    output_fragment_observed = args.output_fragment_observed
    output_fragment_expected = args.output_fragment_expected
    correction_weights_file = args.correction_weights
    fragment_correction_weights_file = args.fragment_correction_weights
    correction_method = args.correction_method
    fragment_min = args.fragment_min
    fragment_max = args.fragment_max
    chromosomes = args.chromosomes
    num_cores = args.num_cores
    twobit_file = args.twobit
    num_simulations = args.num_simulations
    blacklist_file = args.blacklist
    output_bam_file = args.output_bam

    tmp_dir = "tmp_bam_dir"
    os.makedirs(tmp_dir, exist_ok=True)

    logging.info(f"Loading blacklist file: {blacklist_file} ...")
    blacklist_trees = load_blacklist(blacklist_file)
    logging.info(f"Blacklist loaded ...")
    
    logging.info(f"Precomputing allowed intervals for chromosomes ...")
    allowed_intervals_dict = precompute_allowed_intervals(chromosomes, twobit_file, blacklist_trees)
    logging.info(f"Allowed intervals precomputed ...")

    logging.info(f"Processing file: {bam_file} ...")
    bam_results, bam_fragment_counts = parallel_process_bam(bam_file, chromosomes, fragment_min, fragment_max, num_cores)
    logging.info(f"BAM processing complete ...")
    
    observed_motifs = {frag_len: list(motifs.keys()) for frag_len, motifs in bam_results.items()}
    
    logging.info(f"Simulating {num_simulations} end motifs and fragment lengths from genome ...")
    sim_results, sim_fragment_counts = parallel_simulate_motifs(twobit_file, chromosomes, bam_fragment_counts, observed_motifs, num_simulations, num_cores, allowed_intervals_dict)
    logging.info(f"Simulations complete ...")

    if output_observed:
        with open(output_observed, 'w') as f_obs:
            f_obs.write("FragLength\tEndMotif\tCount\n")
            for fragment_length, motifs in sorted(bam_results.items()):
                for motif, count in sorted(motifs.items()):
                    f_obs.write(f"{fragment_length}\t{motif}\t{count}\n")
    
    if output_fragment_observed:
        with open(output_fragment_observed, 'w') as f_frag_obs:
            f_frag_obs.write("FragLength\tCount\n")
            for fragment_length, count in sorted(bam_fragment_counts.items()):
                f_frag_obs.write(f"{fragment_length}\t{count}\n")
    
    if output_expected:
        with open(output_expected, 'w') as f_exp:
            f_exp.write("FragLength\tEndMotif\tCount\n")
            for fragment_length, motifs in sorted(sim_results.items()):
                for motif, count in sorted(motifs.items()):
                    f_exp.write(f"{fragment_length}\t{motif}\t{count}\n")
    
    if output_fragment_expected:
        with open(output_fragment_expected, 'w') as f_frag_exp:
            f_frag_exp.write("FragLength\tCount\n")
            for fragment_length, count in sorted(sim_fragment_counts.items()):
                f_frag_exp.write(f"{fragment_length}\t{count}\n")
    
    if correction_method == 'simple':
        logging.info("Calculating simple correction weights ...")
        motif_weights = calculate_simple_weights(bam_results, sim_results)
        fragment_weights = calculate_simple_fragment_weights(bam_fragment_counts, sim_fragment_counts)
    elif correction_method == 'glm':
        logging.info("Calculating GLM correction weights ...")
        motif_weights = calculate_glm_weights(bam_results, sim_results)
        fragment_weights = calculate_glm_fragment_weights(bam_fragment_counts, sim_fragment_counts)
    
    logging.info("Applying extreme weight limitation ...")
    motif_weights = limit_extreme_weights(motif_weights, args.stringency)
    fragment_weights = limit_extreme_fragment_weights(fragment_weights, args.stringency)

    logging.info("Applying Gaussian smoothing ...")
    motif_weights = smooth_weights(motif_weights, args.smoothing_sigma)
    fragment_weights = smooth_fragment_weights(fragment_weights, args.smoothing_sigma)

    logging.info("Writing correction weights to files ...")
    write_weights_to_tsv(motif_weights, correction_weights_file)
    write_fragment_weights_to_tsv(fragment_weights, fragment_correction_weights_file)
    logging.info("Correction weights written to files.")

    if output_bam_file:
        logging.info("Adding tags to BAM file ...")
        parallel_process_bam_with_tags(bam_file, output_bam_file, chromosomes, motif_weights, fragment_weights, fragment_min, fragment_max, num_cores, tmp_dir)
        logging.info("Tagged BAM file written.")

    logging.info("Complete ...")
    shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    main()