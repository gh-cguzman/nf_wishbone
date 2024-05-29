#!/usr/bin/env python
import sys
import pysam
from collections import defaultdict
import itertools
import argparse

# Parameters:
insert_size_range = (110, 230)

def filter_read(read):
    """Filters reads based on mapping quality, CIGAR string, and reference name."""
    return (
        read.mapping_quality >= 60
        and read.cigarstring == "101M"
        and read.reference_name in ("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")
    )

def process_bam(inBam, EM_out, FEMS_out, QC_out):
    """Processes a GCparagon corrected BAM file and generates the required output files."""
    sample_id_tmp = inBam.split("/")[-1]
    sample_id = f"{sample_id_tmp.split('.')[0]}_GC"
    print(sample_id)

    bamfile = pysam.AlignmentFile(inBam, "rb")

    EM_count = defaultdict(int)
    FEMS_count = defaultdict(int)

    rc = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}

    processed_reads_counter = 0
    bam_filtered_reads_counter = 0
    insert_size_filtered_counter = 0
    bad_flag_filtered_counter = 0
    N_filtered_reads_counter = 0

    for read in bamfile:
        processed_reads_counter += 1

        if not filter_read(read):
            bam_filtered_reads_counter += 1
            continue

        insert_size = abs(read.template_length)
        if not (insert_size_range[0] <= insert_size <= insert_size_range[1]):
            insert_size_filtered_counter += 1
            continue

        seq = read.seq.upper()

        if read.flag in (99, 163):
            end_motif = seq[:4]
        elif read.flag in (83, 147):
            end_motif = ''.join(rc[base] for base in reversed(seq[-4:]))
        else:
            bad_flag_filtered_counter += 1
            continue

        if "N" in end_motif:
            N_filtered_reads_counter += 1
            continue

        gc_tag = read.get_tag("GC") if read.has_tag("GC") else 1

        key = (end_motif, insert_size)
        EM_count[end_motif] += gc_tag
        FEMS_count[key] += gc_tag

    bamfile.close()

    E = ["".join(p) for p in itertools.product('ACGT', repeat=4)]
    I = list(range(110, 231))

    with open(EM_out, 'wt') as f:
        f.write("sample_ID\t" + '\t'.join(E) + '\n')
    with open(FEMS_out, 'wt') as f:
        f.write("sample_ID\t" + '\t'.join(f"{e}-{i}" for e in E for i in I) + '\n')

    for count, outfile in [(EM_count, EM_out), (FEMS_count, FEMS_out)]:
        with open(outfile, 'a') as f:
            if outfile == EM_out:
                f.write(sample_id + '\t' + '\t'.join(str(count[e]) for e in E) + '\n')
            if outfile == FEMS_out:
                f.write(sample_id + '\t' + '\t'.join(str(count[(e, i)]) for e in E for i in I) + '\n')

    with open(QC_out, 'wt') as f:
        f.write("sample_ID\tprocessed_reads\tbam_filtered_reads\tinsert_size_filtered_reads\tbad_flag_filtered_reads\tN_filtered_reads\tPassing_FEMS_reads\n")
        f.write(
            f"{sample_id}\t{processed_reads_counter}\t{bam_filtered_reads_counter}\t{insert_size_filtered_counter}\t{bad_flag_filtered_counter}\t{N_filtered_reads_counter}\t{sum(FEMS_count.values())}\n"
        )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process BAM file and generate EndMotif and FEMS counts.')
    parser.add_argument('--bam', type=str, required=True, help='Path to the BAM file')
    parser.add_argument('--out_em', type=str, required=True, help='Path to the EndMotif output file')
    parser.add_argument('--out_fems', type=str, required=True, help='Path to the FEMS output file')
    parser.add_argument('--out_qc', type=str, required=True, help='Path to the QC output file')
    args = parser.parse_args()

    print(f"\nCarrying out FEMS processing ...\n")
    process_bam(args.bam, args.out_em, args.out_fems, args.out_qc)
