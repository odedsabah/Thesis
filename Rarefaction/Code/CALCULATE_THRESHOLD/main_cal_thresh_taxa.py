#!/usr/bin/python3
'''
Usage example:
python3 main_cal_thresh_taxa.py /data2/Human/HMP.phase1/Stool/raw.d/SAMN00040286/SRS019068.denovo_duplicates_marked.trimmed.1.fastq.gz /data2/Human/HMP.phase1/Stool/raw.d/Stool.read-stats.txt /data1/Oded.Sabah/metaanalysis/SAMN00040286
'''
import sys
import os
import argparse
from fastq_selector import FastqProcessor
from fastq_splitter import FastqSplitter
from Stat_from_select import MetaphlanAnalysis
from Stat_from_split import MetaphlanAnalysis_from_split
from concurrent.futures import ProcessPoolExecutor

def parse_arguments():
    parser = argparse.ArgumentParser(description='Calculate threshold taxa from FASTQ files')
    parser.add_argument('fastq_file', help='Path to the FASTQ file')
    parser.add_argument('read_stats', help='Path to the Read Stats file')
    parser.add_argument('output_dir', help='Directory name for output')
    parser.add_argument('--thresholds', nargs='+', type=float, default=[0, 0.001, 0.01, 0.1, 1],
                        help='Thresholds for analysis (default: [0, 0.001, 0.01, 0.1, 1])')
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    # Create output directory if it doesn't exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    with ProcessPoolExecutor() as executor:
        processor = FastqProcessor(fastq_file=args.fastq_file, read_stats=args.read_stats,
                                directory_name=args.output_dir)
        splitter = FastqSplitter(fastq_file=args.fastq_file, read_stats=args.read_stats,
                              directory_name=args.output_dir)

        processor_future = executor.submit(processor.process)
        splitter_future = executor.submit(splitter.split)

        # Wait for both processes to finish and get their results
        path_selector_out = processor_future.result()
        path_splitter_out = splitter_future.result()

        analysis_future1 = executor.submit(MetaphlanAnalysis, path_mp4_files=path_selector_out,
                                        directory_name=args.output_dir)
        analysis_future2 = executor.submit(MetaphlanAnalysis_from_split, path_mp4_files=path_splitter_out,
                                        directory_name=args.output_dir)

        analysis1 = analysis_future1.result()
        analysis1.run_analysis(thresholds=args.thresholds)

        analysis2 = analysis_future2.result()
        analysis2.run_analysis_from_split(thresholds=args.thresholds)

if __name__ == '__main__':
    main()
