#!/usr/bin/python3

import math
import pandas as pd
import gzip
import os
import logging

'''
This module splits a FASTQ file into multiple groups and runs Metaphlan on each group.

The script takes three required inputs:
1. fastq_file: The path to the FASTQ file.
2. read_stats: The path to the read stats file.
3. directory_name: The directory where results will be stored.
'''

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class FastqSplitter:
    def __init__(self, fastq_file, read_stats, directory_name):
        self.fastq_file = fastq_file
        self.read_stats = read_stats
        self.directory_name = directory_name
        self.sample_size, self.id_file = self.cal_sample_size()
        self.fout_list = []
        self.list2mp4 = []
        self.output_dir = f"{self.directory_name}/mp4.split.{self.id_file}"
        if not os.path.isdir(self.output_dir):
            os.mkdir(self.output_dir)
            logging.info(f"Created output directory: {self.output_dir}")
        logging.info(f"Initialized FastqSplitter for {fastq_file}")

    def prepare_fastq_file(self, jump):
        try:
            with gzip.open(self.fastq_file, "rt") as fin:
                for i, line in enumerate(fin):
                    line += fin.readline()
                    line += fin.readline()
                    line += fin.readline()
                    group = i % 10
                    if (i - group) % jump == 0:
                        self.fout_list[group].write(line)
                        if i == self.sample_size:
                            break
            logging.info(f"Split FASTQ file with jump={jump}")
        except Exception as e:
            logging.error(f"Error preparing split FASTQ file: {e}")
            raise

    def cal_sample_size(self):
        try:
            id_file = self.fastq_file.split('/')[-2]
            df_size = pd.read_csv(self.read_stats, delimiter='\t', skiprows=2)
            size_by_id = df_size[df_size.iloc[:, 0].str.startswith(id_file)]
            if size_by_id.empty:
                raise ValueError(f"No entries found for ID {id_file} in read stats file")
                
            sample_size = int(size_by_id.iloc[:, 3].iloc[0])
            sample_size = math.floor(sample_size / 1e7) * 1e7
            logging.info(f"Calculated sample size: {sample_size} for {id_file}")
            return sample_size, id_file
        except Exception as e:
            logging.error(f"Error calculating sample size: {e}")
            raise

    def split(self):
        logging.info(f"Starting to split {self.fastq_file}")
        num_reads = 0
        while True:
            max_groups = 10
            num_reads += 1e6
            output_files = [f"{self.output_dir}/{self.id_file}_group_{i+1}-{int(num_reads//1e6)}M.fastq.gz" 
                           for i in range(0, max_groups)]

            # Define the number of reads to select and the jump value
            jump = int(self.sample_size / num_reads)
            if jump < max_groups:
                logging.info(f"Reached minimum jump ({jump}), stopping at {num_reads} reads")
                break

            logging.info(f"Processing {num_reads} reads with jump={jump}")
            # Open all output files at the beginning
            self.fout_list = [gzip.open(output_file, "wt") for output_file in output_files]
            self.list2mp4.extend(self.fout_list)
            
            # Split the fastq file into 10 groups
            self.prepare_fastq_file(jump)

            # Close all output files
            for fout in self.fout_list:
                fout.close()

        for file_object in self.list2mp4:
            path_file2mp4 = file_object.name.split(".fastq.gz")[0]
            metaphlan_file_in = f"{file_object.name}"
            metaphlan_file_out = f"{path_file2mp4}.txt"
            cmd = f'/data1/software/metaphlan/run-metaphlan.sh {metaphlan_file_in} {metaphlan_file_out} 40 > {metaphlan_file_out}.stdout'
            logging.info(f"Running Metaphlan: {cmd}")
            os.system(cmd)
            
        return self.output_dir
