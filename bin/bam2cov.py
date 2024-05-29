#!/usr/bin/env python
import sys
import pysam
from collections import defaultdict
import argparse
from intervaltree import Interval, IntervalTree
import re
import multiprocessing as mp

def load_blacklist(blacklist_file):
    """Loads the blacklist BED file and returns an IntervalTree for each chromosome."""
    blacklist = defaultdict(IntervalTree)
    with open(blacklist_file, 'r') as f:
        for line in f:
            if line.strip():
                chrom, start, end = line.strip().split()[:3]
                blacklist[chrom].add(Interval(int(start), int(end)))
    return blacklist

def load_regions(region_file):
    """Loads the region BED file and returns an IntervalTree for each chromosome."""
    region = defaultdict(IntervalTree)
    with open(region_file, 'r') as f:
        for line in f:
            if line.strip():
                chrom, start, end = line.strip().split()[:3]
                region[chrom].add(Interval(int(start), int(end)))
    return region

def filter_reads(read, blacklist, region, min_mapq=60, min_length=110, max_length=230):
    if read.mapping_quality < min_mapq:
        return False
    if read.flag not in {99, 147, 83, 163}:
        return False
    if read.reference_name not in region:
        return False
    if not min_length <= abs(read.template_length) <= max_length:
        return False
    if read.reference_name in blacklist and blacklist[read.reference_name].overlaps(read.reference_start, read.reference_end):
        return False
    return True

def chrom_sort_key(chrom):
    """Sort key function to handle chromosome names lexicographically and numerically."""
    match = re.match(r'chr(\d+|[XYM])$', chrom)
    if match:
        value = match.group(1)
        if value.isdigit():
            return (int(value),)
        else:
            return (value,)
    return (chrom,)

def worker(bam_file, regions, blacklist, chrom, return_dict, idx):
    """Worker function to process a specific chromosome in the BAM file."""
    read_counts = defaultdict(int)
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(chrom):
            if filter_reads(read, blacklist, regions):
                for region_interval in regions[read.reference_name].overlap(read.reference_start, read.reference_end):
                    region_start, region_end = region_interval.begin, region_interval.end
                    region_key = (read.reference_name, region_start, region_end)
                    gc_tag = read.get_tag("GC") if read.has_tag("GC") else 1
                    read_counts[region_key] += gc_tag
    return_dict[idx] = read_counts

def process_bam_parallel(bam_file, regions, blacklist, outcounts, outlog, threads, sample):
    """Process the BAM file in parallel and count overlaps with regions."""
    read_counts = defaultdict(int)

    manager = mp.Manager()
    return_dict = manager.dict()
    jobs = []

    chromosomes = list(regions.keys())
    num_chroms = len(chromosomes)
    chunk_size = max(1, num_chroms // threads)

    for i in range(0, num_chroms, chunk_size):
        chunk_chroms = chromosomes[i:i + chunk_size]
        for chrom in chunk_chroms:
            p = mp.Process(target=worker, args=(bam_file, regions, blacklist, chrom, return_dict, chrom))
            jobs.append(p)
            p.start()

    for proc in jobs:
        proc.join()

    # Combine results from all processes
    for result in return_dict.values():
        for key, count in result.items():
            read_counts[key] += count

    # Ensure that all regions are accounted for even if they have zero counts
    all_regions = {(chrom, interval.begin, interval.end) for chrom, intervals in regions.items() for interval in intervals}
    for region in all_regions:
        if region not in read_counts:
            read_counts[region] = 0

    with open(outcounts, "w") as out_counts, open(outlog, "w") as out_log:
        out_log.write(f"Sample: {sample}\n")
        for (chrom, start, end), count in sorted(read_counts.items(), key=lambda x: (chrom_sort_key(x[0][0]), x[0][1], x[0][2])):
            out_counts.write(f"{chrom}\t{start}\t{end}\t{count}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process BAM file and count overlaps.")
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--regions", required=True, help="Regions BED file")
    parser.add_argument("--blacklist", required=True, help="Blacklist BED file")
    parser.add_argument("--outcounts", required=True, help="Output file for read count summary")
    parser.add_argument("--outlog", required=True, help="Output log file")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads for multithreading")
    parser.add_argument("--sample", required=True, help="Sample name")
    args = parser.parse_args()

    print(f"\nRunning sample {args.sample} ...")

    regions = load_regions(args.regions)
    blacklist = load_blacklist(args.blacklist)
    process_bam_parallel(args.bam, regions, blacklist, args.outcounts, args.outlog, args.threads, args.sample)
