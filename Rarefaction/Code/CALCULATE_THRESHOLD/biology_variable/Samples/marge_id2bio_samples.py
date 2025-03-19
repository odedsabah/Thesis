#!/usr/bin/python3
import pandas as pd
import re
import sys
import subprocess
import os
import time


"""
This script is designed to process and cross-reference data from two distinct data sources related to the study published in https://doi.org/10.1038/nbt.3960.

Process Overview:
1. Input:
   - The script takes two DataFrames (df) as input.
   - The first DataFrame is a summary of samples obtained from the European Nucleotide Archive (ENA) at https://www.ebi.ac.uk/ena/browser/view/PRJEB14847. This DataFrame contains detailed information about various samples.
   - The second DataFrame is sourced from the supplementary materials of the mentioned article. This DataFrame provides additional context and description for each sample.

2. Cross-referencing:
   - The script cross-references these two DataFrames based on the sample IDs. This step involves matching the sample IDs from the ENA summary DataFrame with those in the article's supplementary DataFrame to integrate the data.

3. Output:
   - The output of the script is a consolidated DataFrame. This DataFrame combines the matched data from both sources, providing a comprehensive view of each sample, including its detailed description and contextual information from the article's supplementary material.

Purpose:
- The script aims to filter and download 41 unique samples based on maximum base count from a total of 21 protocols, each having two stool samples (A & B). It addresses the scenario where the total count is 41 instead of 42 due to a mismatch in one of the samples.
"""

def concatenate_dataframes(df1, df2):
    # Merge the DataFrames based on "Sample" and "sample_alias" columns
    merged_df = df1.merge(df2, left_on="Sample", right_on="sample_alias", how="inner")
    # Create a new column 'sample_gro' containing the first letter of 'Sample'
    merged_df['sample_gro'] = merged_df['Sample'].apply(lambda x: x[0])

    # Group by "Lab" and "sample_gro" and get the maximum value in "base_count" for each group
    def custom_group(group):
        return group[group["base_count"] == group["base_count"].max()]

    grouped_df = merged_df.groupby([ "sample_gro","Lab"]).apply(custom_group)

    return grouped_df.reset_index(drop=True)

def replace_last_underscore(alias):
    parts = alias.rsplit("_", maxsplit=1)
    if len(parts) > 1:
        return parts[0] + "-" + parts[1]
    else:
        return alias

def main(path_id_files, path_filereport_read):
    csv_data = pd.read_csv(path_id_files)
    txt_data = pd.read_csv(path_filereport_read, sep='\t')

    # Extract relevant parts and find matches
    csv_data['Sample'] = csv_data['Sample'].apply(lambda x: '_'.join(x.split('_')[:2]))

    # Clean 'sample_alias' in txt_data
    txt_data['sample_alias'] = txt_data['sample_alias'].apply(lambda x: replace_last_underscore(x))
    txt_data['sample_alias'] = txt_data['sample_alias'].apply(
        lambda x: x.split("_")[-1].replace(" ", "_").replace("-", "_").split("749906_")[-1])
    concatenate = concatenate_dataframes(csv_data, txt_data)
    concatenate = concatenate.drop_duplicates()
    # concatenate.to_csv("/home/odeds/metaanalysis/Rarefaction/Code/CALCULATE_THRESHOLD/biology_variable/Samples/read-sample.csv")
    BASE_DEST_FOLDER = "/data1/Oded.Sabah/metaanalysis/biology_variable/raw.d"

    # Iterate over the DataFrame and download files
    for index, row in concatenate.iterrows():
        if pd.notna(row['fastq_ftp']) and pd.notna(row['sample_accession']):
            # Create a directory for this sample
            sample_dir = os.path.join(BASE_DEST_FOLDER, row['sample_accession'])
            os.makedirs(sample_dir, exist_ok=True)

            # Split the FTP URLs
            ftp_urls = row['fastq_ftp'].split(';')

            # Download each file using wget
            for url in ftp_urls:
                subprocess.run(['wget', '-P', sample_dir, url])
                print(f"Downloaded {url} to {sample_dir}")
                time.sleep(2)  # 2-second delay between downloads

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("\nUsage: " + sys.argv[0] + " <Path to file ID csv> <Path to filereport read txt>\n\n")
    else:
        main(sys.argv[1], sys.argv[2])