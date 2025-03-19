#!/usr/bin/python3
#exemple
'''
python3 main_cal_thresh_taxa.py /data2/Human/HMP.phase1/Stool/raw.d/SAMN00040286/SRS019068.denovo_duplicates_marked.trimmed.1.fastq.gz /data2/Human/HMP.phase1/Stool/raw.d/Stool.read-stats.txt /data1/Oded.Sabah/metaanalysis/SAMN00040286
'''
import sys
from fastq_selector import FastqProcessor
from fastq_splitter import FastqSplitter
from Stat_from_select import MetaphlanAnalysis
from Stat_from_split import MetaphlanAnalysis_from_split
from concurrent.futures import ProcessPoolExecutor, wait

def main():
    if len(sys.argv) != 4:
        print("\nUsage: " + sys.argv[0] + " <Path to the FASTQ file> <Path to the Read Stats file> <Directory name out> \n\n")

    Path_to_the_FASTQ_file = sys.argv[1]
    Path_to_the_Read_Stats_file = sys.argv[2]
    directory_name = sys.argv[3]

    with ProcessPoolExecutor() as executor:
        processor = FastqProcessor(fastq_file=Path_to_the_FASTQ_file, read_stats=Path_to_the_Read_Stats_file,
                                   directory_name=directory_name)
        splitter = FastqSplitter(fastq_file=Path_to_the_FASTQ_file, read_stats=Path_to_the_Read_Stats_file,
                                 directory_name=directory_name)

        processor_future = executor.submit(processor.process)
        splitter_future = executor.submit(splitter.split)

        # Wait for both processes to finish and get their results
        path_selector_out = processor_future.result()
        path_splitter_out = splitter_future.result()


        # path_selector_out = '/data1/Oded.Sabah/metaanalysis/Rarefaction_out/SAMN00143093/mp4.depth.SAMN00143093'
        # path_splitter_out = '/data1/Oded.Sabah/metaanalysis/Rarefaction_out/SAMN00143093/mp4.split.SAMN00143093'
        # directory_name = '/data1/Oded.Sabah/metaanalysis/Rarefaction_out/SAMN00143093'
        #
        threshold = [0, 0.001, 0.01, 0.1, 1]
        analysis_future1 = executor.submit(MetaphlanAnalysis, path_mp4_files=path_selector_out,
                                           directory_name=directory_name)
        analysis_future2 = executor.submit(MetaphlanAnalysis_from_split, path_mp4_files=path_splitter_out,
                                           directory_name=directory_name)

        analysis1 = analysis_future1.result()
        analysis1.run_analysis(thresholds=threshold)

        analysis2 = analysis_future2.result()
        analysis2.run_analysis_from_split(thresholds=threshold)

if __name__ == '__main__':
    main()
