#!/usr/bin/python3

# Upload packages
import gzip
import os
import pandas as pd
import math

'''
In the code, `id_file` is derived from the input FASTQ file's ID, and `num_reads` represents the number of reads selected for each iteration.
The output files are created for each iteration of the loop in the `main` function, 
where subsets of reads are processed and analyzed by Metaphlan.

The script takes two command-line arguments as input:
1. path_fastq_file: The path to the FASTQ file.
2. read_stats: The path to the read stats file (with end ".txt").
'''

class FastqProcessor:
    def __init__(self, fastq_file, read_stats, directory_name):
        self.fastq_file = fastq_file
        self.read_stats = read_stats
        self.sample_size, self.id_file = self.cal_sample_size()
        self.dict_selected_read = self.selected_reads()
        self.directory_name = directory_name

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

    def cal_sample_size(self):
        id_file = self.fastq_file.split('/')[-2]
        df_size = pd.read_csv(self.read_stats, delimiter='\t', skiprows=1)
        size_by_id = df_size[df_size.iloc[:, 0].str.startswith(id_file)]
        sample_size = int(size_by_id.iloc[:, 3].iloc[0])
        sample_size = math.floor(sample_size / 1e7) * 1e7
        return sample_size, id_file

    def process(self):
        global path2process
        id_file = f"{self.directory_name}/mp4.depth.{self.fastq_file.split('/')[-2]}/" + self.fastq_file.split('/')[-2]
        if not os.path.isdir(f"{self.directory_name}/mp4.depth.{self.fastq_file.split('/')[-2]}"):
            os.mkdir(f"{self.directory_name}/mp4.depth.{self.fastq_file.split('/')[-2]}")
        for i, (num_reads, jump) in enumerate(self.dict_selected_read.items()):
            fastq_file  = f"{id_file}.{int(num_reads // 1e6)}M.fastq.gz"
            fastq_file_out = f"{id_file}.{int(num_reads // 1e6)}M.fastq.gz"
            metaphlan_file_out = f"{id_file}.{int(num_reads // 1e6)}M.txt"
            path2process = (f"{self.directory_name}/mp4.depth.{self.fastq_file.split('/')[-2]}")
            self.prepare_fastq_file(fastq_file, jump, num_reads)
            cmd = (f'/data1/software/metaphlan/run-metaphlan.sh {fastq_file_out} {metaphlan_file_out}'
                      f' 40 > {metaphlan_file_out}.stdout')
            os.system(cmd)
            # os.system(f"rm -f {fastq_file_out}")
        return path2process


