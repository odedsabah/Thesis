# #!/usr/bin/python3

# '''
# python3 ~/metaanalysis/Rarefaction/Code/CALCULATE_THRESHOLD/biology_variable/Shuffle_thresh2biology_variable.py /data1/Human/MetaHit/raw.d/ /data1/Oded.Sabah/metaanalysis/biology_variable/
# '''

import gzip
import os
import random
import sys
from Bio import SeqIO

class RemoveReads:
    def __init__(self, base_path, directory_name):
        self.base_path = os.path.expanduser(base_path)
        self.directory_name = directory_name
        self.allowed_ids = ["SAMEA728859", "SAMEA728763", "SAMEA728929",
            "SAMEA728734", "SAMEA728680", "SAMEA728611",
            "SAMEA728862", "SAMEA728916"] #MetaHit samples
        n_serial = (list(range(1, 10)) + list(range(10, 19, 2)) +
                    list(range(20, 49, 5)) + list(range(50, 99, 10)) +
                    list(range(100, 201, 100)))
        self.sample_size = [i * 1e6 for i in n_serial]

    def find_fastq_file(self, sample_id):
        sample_path = os.path.join(self.base_path, sample_id)
        if os.path.isdir(sample_path):
            for file_name in sorted(os.listdir(sample_path)):  # Sort to get consistent results
                if file_name.endswith('.fastq') or file_name.endswith('.fastq.gz'):
                    return os.path.join(sample_path, file_name)
        return None

    def randomly_subsample_fastq(self, sample_id, fastq_file):
        if fastq_file.endswith('.gz'):
            open_func = gzip.open
            mode = 'rt'  # Read as text mode
        else:
            open_func = open
            mode = 'r'

        with open_func(fastq_file, mode) as handle:
            records = list(SeqIO.parse(handle, "fastq"))
            total_reads = len(records)

        for sample_size in self.sample_size:
            if sample_size > total_reads:
                print(
                    f"Skipping subsampling for {sample_size} as it exceeds total reads ({total_reads}) in {fastq_file}")
                continue

            output_filename = f"{sample_id}.subsampled_{int(sample_size // 1e6)}M.fastq.gz"
            output_fastq = os.path.join(self.directory_name, sample_id, output_filename)

            with gzip.open(output_fastq, 'wt') as output_handle:
                for record in records[:int(sample_size)]:
                    SeqIO.write(record, output_handle, "fastq")

            yield output_fastq

    def process(self):
        for sample_id in self.allowed_ids:
            fastq_file = self.find_fastq_file(sample_id)
            if fastq_file:
                print(f"Processing FASTQ file: {fastq_file}")
                id_file_dir = os.path.join(self.directory_name, sample_id)
                if not os.path.isdir(id_file_dir):
                    os.makedirs(id_file_dir)

                for fastq_file_out in self.randomly_subsample_fastq(sample_id, fastq_file):
                    metaphlan_file_out = os.path.join(id_file_dir, f'{os.path.basename(fastq_file_out.split(".fastq.gz")[0])}.txt')

                    cmd = (f'/data1/software/metaphlan/run-metaphlan.sh {fastq_file_out} {metaphlan_file_out}'
                           f' 40 > {metaphlan_file_out}.stdout')
                    os.system(cmd)
            else:
                print(f"No FASTQ file found for {sample_id}")

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <Base path to the FASTQ files> <Directory name out>")
        sys.exit(1)

    base_path = sys.argv[1]
    directory_name_out = sys.argv[2]

    remover = RemoveReads(base_path, directory_name_out)
    remover.process()

if __name__ == '__main__':
    main()
