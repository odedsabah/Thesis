#!/usr/bin/python3

import sys
import os


def id_files(path_fastq_file):
    path_id = {}
    # Iterate over directories in the given path
    sample_dirs = [d for d in os.listdir(path_fastq_file) if os.path.isdir(os.path.join(path_fastq_file, d))]

    for dir_name in sample_dirs:
        sample_dir_path = os.path.join(path_fastq_file, dir_name)
        # Iterate over files in each sample directory
        samples = [f for f in os.listdir(sample_dir_path) if f.endswith(".gz")]

        for sample in samples:
            id_sample = sample.split(".fastq.gz")[0]
            path_sample = os.path.join(sample_dir_path, sample)
            path_id[id_sample] = path_sample

    return path_id


def main():
    if len(sys.argv) != 2:
        quit("\nUsage: " + sys.argv[0] + " <Path_MP4_files> \n\n")

    path_MP4_files = sys.argv[1]

    path_id = id_files(path_MP4_files)

    for id_file, path_file in path_id.items():
        metaphlan_file_out = f"{path_file.split('.fastq.gz')[0]}.txt"
        print(metaphlan_file_out)
        cmd = f'/data1/software/metaphlan/run-metaphlan.sh {path_file} {metaphlan_file_out} 40 > {metaphlan_file_out}.stdout'
        os.system(cmd)

if __name__ == '__main__':
    main()
