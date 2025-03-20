#!/usr/bin/python3

# Upload packages
import gzip
import os
import pandas as pd
import math
import logging

'''
This module selects subsets of reads from a FASTQ file at different depths and runs Metaphlan on them.

The script takes two required inputs:
1. path_fastq_file: The path to the FASTQ file.
2. read_stats: The path to the read stats file (with end ".txt").
3. directory_name: The directory where results will be stored.
'''

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class FastqProcessor:
    def __init__(self, fastq_file, read_stats, directory_name):
        self.fastq_file = fastq_file
        self.read_stats = read_stats
        self.directory_name = directory_name
        self.sample_size, self.id_file = self.cal_sample_size()
        self.dict_selected_read = self.selected_reads()
        logging.info(f"Initialized FastqProcessor for {fastq_file}")

    def selected_reads(self):
        n_serial = (list(range(1, 10)) + list(range(10, 19, 2)) +
                    list(range(20, 49, 5)) + list(range(50, 99, 10)) +
                    list(range(100, 1001, 100)))
        n_serial = [i * 1e6 for i in n_serial]
        selected_read = [self.sample_size // (i if i != 0 else 1) for i in n_serial]
        dict_selected_read = dict(zip(n_serial, selected_read))
        dict_selected_read = {k: v for k, v in dict_selected_read.items() if v != 0}
        return dict_selected_read

    def prepare_fastq_file(self, fastq_file_out, jump, num_reads):
        try:
            with gzip.open(self.fastq_file, "rt") as fin, gzip.open(fastq_file_out, "wt") as fout:
                counter = 0
                for line in fin:
                    line += fin.readline()
                    line += fin.readline()
                    line += fin.readline()
                    counter += 1
                    if counter % jump == 0:
                        fout.write(line)
                        if counter // jump == num_reads:
                            break
            logging.info(f"Created {fastq_file_out} with {num_reads} reads (jump={jump})")
        except Exception as e:
            logging.error(f"Error preparing FASTQ file: {e}")
            raise

    def cal_sample_size(self):
        try:
            id_file = self.fastq_file.split('/')[-2]
            df_size = pd.read_csv(self.read_stats, delimiter='\t', skiprows=1)
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

    def process(self):
        output_dir = f"{self.directory_name}/mp4.depth.{self.id_file}"
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
            logging.info(f"Created output directory: {output_dir}")
            
        id_file = f"{output_dir}/{self.id_file}"
        
        for num_reads, jump in self.dict_selected_read.items():
            fastq_file_out = f"{id_file}.{int(num_reads // 1e6)}M.fastq.gz"
            metaphlan_file_out = f"{id_file}.{int(num_reads // 1e6)}M.txt"
            
            logging.info(f"Processing {num_reads} reads...")
            self.prepare_fastq_file(fastq_file_out, jump, num_reads)
            
            cmd = (f'/data1/software/metaphlan/run-metaphlan.sh {fastq_file_out} {metaphlan_file_out}'
                  f' 40 > {metaphlan_file_out}.stdout')
            logging.info(f"Running Metaphlan: {cmd}")
            os.system(cmd)
            
            # Uncomment to remove the fastq file after processing
            # os.system(f"rm -f {fastq_file_out}")
            
        return output_dir


