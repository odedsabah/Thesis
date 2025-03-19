# #!/usr/bin/python3

import math
import pandas as pd
import gzip
import sys
import os

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

    def prepare_fastq_file(self, jump):
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

    def cal_sample_size(self):
        id_file = self.fastq_file.split('/')[-2]
        df_size = pd.read_csv(self.read_stats, delimiter='\t', skiprows=2)
        size_by_id = df_size[df_size.iloc[:, 0].str.startswith(id_file)]
        sample_size = int(size_by_id.iloc[:, 3].iloc[0])
        sample_size = math.floor(sample_size / 1e7) * 1e7
        return sample_size, id_file

    def split(self):
        global path2process
        path2process = f"{self.directory_name}/mp4.split.{self.id_file}"
        num_reads = 0
        while True:
            max_groups = 10
            num_reads += 1e6
            output_files = [f"{self.output_dir}/{self.id_file}_group_{i+1}-{int(num_reads//1e6)}M.fastq.gz" for i in range(0, max_groups)]

            # Define the number of reads to select and the jump value
            jump = int(self.sample_size / num_reads)
            if jump < max_groups:
                break

            # Open all output files at the beginning
            self.fout_list = [gzip.open(output_file, "wt") for output_file in output_files]
            self.list2mp4.extend(self.fout_list)
            # Split the fastq file into 10 groups of 1 million reads each
            # Select the sequences for the group and write them to the output files
            self.prepare_fastq_file(jump)

            #Close all output files
            for fout in self.fout_list:
                fout.close()

        for file_object in self.list2mp4:
            path_file2mp4 = file_object.name.split(".fastq.gz")[0]
            metaphlan_file_in = f"{file_object.name}"
            metaphlan_file_out = f"{path_file2mp4}.txt"
            cmd = f'/data1/software/metaphlan/run-metaphlan.sh {metaphlan_file_in} {metaphlan_file_out} 40 > {metaphlan_file_out}.stdout'
            os.system(cmd)
        return path2process
