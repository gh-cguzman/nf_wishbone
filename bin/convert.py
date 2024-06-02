#!/usr/bin/env python
import pandas as pd
import argparse

def read_cov(df_cov_path, sep='\t', header=None, index_col=None):
    df = pd.read_csv(df_cov_path, sep=sep, header=header, index_col=index_col)
    return df

def convert(df_cov, sample_id, cov_loci_cols=[0, 1, 2], cov_value_col=3):
    # Construct loci from the specified columns of df_cov
    loci = df_cov[cov_loci_cols[0]].astype(str) + "_" + df_cov[cov_loci_cols[1]].astype(str) + "_" + df_cov[cov_loci_cols[2]].astype(str)
    final_df = pd.DataFrame([df_cov[cov_value_col].values], columns=loci)
    
    # Set the index of final_df
    final_df.index = [sample_id]

    return final_df

def main():
    parser = argparse.ArgumentParser(description='Process some files.')
    parser.add_argument('--cov', type=str, required=True, help='Path to the cov file')
    parser.add_argument('--output', type=str, required=True, help='Path to the output TSV file')
    parser.add_argument('--sample_id', type=str, required=True, help='Sample ID')
    parser.add_argument('--cov_loci_cols', type=int, nargs=3, default=[0, 1, 2], help='Column indices for loci in cov file')
    parser.add_argument('--cov_value_col', type=int, default=3, help='Column index for values in cov file')

    args = parser.parse_args()

    df_cov = read_cov(args.cov)
    print(f'\nRunning {file_name} ...')

    final_df = convert(df_cov, args.sample_id, 
                       args.cov_loci_cols, args.cov_value_col)

    final_df.to_csv(args.output, sep='\t', index=True)

if __name__ == "__main__":
    main()
