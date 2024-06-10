#!/usr/bin/env python
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.fft import fft
from scipy.signal import find_peaks
import logging
import os

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract features from GC-corrected coverage profiles")
    parser.add_argument("--input", required=True, help="Input matrix file from the previous script")
    parser.add_argument("--output", required=True, help="Output TSV file for the extracted features")
    parser.add_argument("--plot_dir", required=True, help="Directory to save the plots")
    return parser.parse_args()

def calculate_features(df, plot_dir):
    positions = np.array([col for col in df.columns if col not in ['Sample_ID', 'BED_File']], dtype=int)
    center_index = np.where(positions == 0)[0][0]
    
    features = []
    for _, row in df.iterrows():
        sample_id = row['Sample_ID']
        bed_file = row['BED_File']
        coverage = row[2:].values.astype(float)

        logging.info(f"Processing sample {sample_id}, bed file {bed_file}")

        # Central coverage
        central_coverage_start = max(center_index - 50, 0)
        central_coverage_end = min(center_index + 51, len(coverage))
        central_coverage = np.mean(coverage[central_coverage_start: central_coverage_end])
        
        # Mean coverage
        mean_coverage_start = max(center_index - 1000, 0)
        mean_coverage_end = min(center_index + 1001, len(coverage))
        mean_coverage = np.mean(coverage[mean_coverage_start: mean_coverage_end])

        # Amplitude (using Fast Fourier Transform)
        fft_start = max(center_index - 960, 0)
        fft_end = min(center_index + 961, len(coverage))
        fft_result = np.abs(fft(coverage[fft_start: fft_end]))
        amplitude = fft_result[10]

        # Global maxima and minima within +/- 250 bp regions
        peaks, _ = find_peaks(coverage, prominence=0.01)
        troughs, _ = find_peaks(-coverage, prominence=0.02)

        # Global maxima and minima right
        global_maxima_right = np.nan
        global_minima_right = np.nan
        for peak in peaks:
            if peak > center_index:
                global_maxima_right = positions[peak]
                break
        for trough in troughs:
            if trough > center_index:
                global_minima_right = positions[trough]
                break

        # Global maxima and minima left
        global_maxima_left = np.nan
        global_minima_left = np.nan
        for peak in reversed(peaks):
            if peak < center_index:
                global_maxima_left = positions[peak]
                break
        for trough in reversed(troughs):
            if trough < center_index:
                global_minima_left = positions[trough]
                break

        # Calculate wavelength
        if len(peaks) > 1:
            wavelengths = np.diff(peaks)
            average_wavelength = np.mean(wavelengths)
        else:
            average_wavelength = 0

        features.append([sample_id, bed_file, central_coverage, mean_coverage, amplitude,
                         global_maxima_right, global_maxima_left,
                         global_minima_right, global_minima_left, average_wavelength])
        
        # Plotting Coverage Profile
        plt.figure(figsize=(5, 5))
        plt.plot(positions, coverage, label="Coverage", linewidth=3, color='lightblue')
        plt.axvline(x=positions[central_coverage_start], color='black', linestyle='--', label="Central Coverage Region Start")
        plt.axvline(x=positions[central_coverage_end-1], color='black', linestyle='--', label="Central Coverage Region End")
        plt.axvline(x=positions[mean_coverage_start], color='darkred', linestyle='--', label="Mean Coverage Region Start")
        plt.axvline(x=positions[mean_coverage_end-1], color='darkred', linestyle='--', label="Mean Coverage Region End")
        if not np.isnan(global_maxima_right):
            plt.plot(global_maxima_right, coverage[positions == global_maxima_right], 'go', label="+1 Nuc (Peak)", markersize=10)
        if not np.isnan(global_maxima_left):
            plt.plot(global_maxima_left, coverage[positions == global_maxima_left], 'bo', label="-1 Nuc (Peak)", markersize=10)
        if not np.isnan(global_minima_right):
            plt.plot(global_minima_right, coverage[positions == global_minima_right], 'ro', label="+1 Nuc (Valley)", markersize=10)
        if not np.isnan(global_minima_left):
            plt.plot(global_minima_left, coverage[positions == global_minima_left], 'mo', label="-1 Nuc (Valley)", markersize=10)
        plt.xlabel("Genomic Position (relative to center)")
        plt.ylabel("Normalized Coverage")
        plt.title(f"Coverage Profile for {sample_id} (Feature: {bed_file})")
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.grid(False)
        os.makedirs(plot_dir, exist_ok=True)
        plt.savefig(f"{plot_dir}/{sample_id}_{bed_file}_coverage_plot.png", bbox_inches='tight', dpi=250)
        plt.close()
        logging.info(f"Saved coverage plot for {sample_id}")

        # Plotting FFT Amplitude
        plt.figure(figsize=(5, 5))
        plt.plot(fft_result, label="FFT Amplitude", linewidth=3, color='lightblue')
        plt.axvline(x=10, color='r', linestyle='--', label="10th Frequency Component")
        plt.xlabel("Frequency")
        plt.ylabel("Amplitude")
        plt.title(f"FFT Amplitude for {sample_id}, {bed_file}")
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.grid(False)
        plt.savefig(f"{plot_dir}/{sample_id}_{bed_file}_fft_plot.png", bbox_inches='tight', dpi=250)
        plt.close()
        logging.info(f"Saved FFT plot for {sample_id}")

    return pd.DataFrame(features, columns=['Sample_ID', 'BED_File', 'Central_Coverage', 'Mean_Coverage', 'Amplitude',
                                           'Global_Maxima_Right', 'Global_Maxima_Left',
                                           'Global_Minima_Right', 'Global_Minima_Left', 'Wavelength'])

def main():
    args = parse_arguments()

    logging.info(f"Reading input matrix from {args.input}")
    df = pd.read_csv(args.input, sep='\t')

    logging.info("Calculating features")
    features_df = calculate_features(df, args.plot_dir)

    logging.info(f"Saving features to {args.output}")
    features_df.to_csv(args.output, sep='\t', index=False)

    logging.info("Feature extraction and plotting completed successfully")

if __name__ == "__main__":
    main()
