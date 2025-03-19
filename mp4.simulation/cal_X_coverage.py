#!/usr/bin/python3

import sys
import os
import subprocess
import gzip
import pandas as pd
from scipy.spatial.distance import braycurtis

'''
The script expects two command-line arguments as inputs:

1. path_fastq_file_X_cov: This is the directory path where .fastq.gz files are located.
 These files contain raw sequencing reads data.
2. Taxa_list_ref: This is a file that contains the reference list of microbial taxa. 
It's expected to be a tab-delimited text file with the taxa names.

The output of this script is a printed DataFrame (stat_species), which contains the following columns:

1. Num_of_species: The number of species in each sample or simulation.
2. Similarity_jac: The Jaccard similarity score for each sample or simulation compared to the reference.
3. Dissimilarity: The Bray-Curtis dissimilarity score for each sample or simulation compared to the reference.
4. l1_score and l2_score: The L1 and L2 scores for each sample or simulation compared to the reference.
5. Min_Abundance and Max_Abundance: The minimum and maximum abundances of each taxon across the data columns. 
'''

def Get_id(fastq_file):
    # list the files in the fastq_file directory and filter for .fastq.gz files
    files_list = [file for file in os.listdir(fastq_file) if file.endswith('.fastq.gz')]
    # create a list of file paths and id_file values using list comprehension
    file_paths = [os.path.join(fastq_file, file_name) for file_name in files_list]
    id_files = [f'{file_name.split("_")[3]}_Coverage' for file_name in files_list]
    # return a tuple with the lists of file paths and id_file values
    return file_paths, id_files

def Create_ref(Taxa_list_ref):
    Taxa_list_ref = pd.read_csv(Taxa_list_ref, header= None, delimiter='\t')
    Taxa_list_ref['Species'] = 's__' + Taxa_list_ref.iloc[:,0].str.split().str[:2].str.join("_")
    Taxa_list_ref['Abundance'] = 10
    Taxa_list_ref = Taxa_list_ref.iloc[:, 1:].set_index('Species').sort_index()
    return Taxa_list_ref

def Mp4_txt_2_csv():
    script_path = '~/MGX_MTX/Predict_diagnosis/MetaPhlan/Sort_Metaphlane.py'
    data_dir = '~/metaanalysis/mp4.simulation/'
    subprocess.run(f'python3 {script_path} {data_dir}', shell=True)
    Taxa_sim = pd.read_csv("~/metaanalysis/species_main_sim.csv")
    Taxa_sim = Taxa_sim[Taxa_sim.iloc[:,0].apply(lambda x: x.startswith('s__'))].set_index('Abundance').sort_index()
    return Taxa_sim

def cal_dis(Taxa_sim,Taxa_list_ref):
    combine2ref = pd.concat([Taxa_sim, Taxa_list_ref], axis= 1).dropna() #Because the same species were not found, the deletion was made.
    combine2ref = combine2ref.rename({'Abundance': 'Ref'}, axis=1)
    # script_path = '~/MGX_MTX/Deep_sequence/Deep_sequence2Metaphlan.py'
    # subprocess.run(f'python3 {script_path} {combine2ref}', shell=True)
    # compute_bray_curtis_dissimilarity(combine2ref,Taxa_stats, ground_truth)
    return combine2ref

def compute_jaccard_scores(df):
    # initialize an empty dictionary to hold the similarity scores
    # df = df.rename({'Abundance':'Ref'}, axis=1)
    similarity_scores = {}
    ground_truth = df.columns[-1]
    # loop over all pairs of columns
    for col in df.columns:
        # compute the Jaccard score between the two columns
        species1 = df.loc[df[col] > 0, :].index
        species2 = df.loc[df[ground_truth] > 0, :].index
        intersection = len(set(species1).intersection(set(species2)))
        union = len(set(species1).union(set(species2)))
        score = intersection / union
        # add the score to the similarity_scores dictionary
        similarity_scores[col] = score
    # create new columns in the dataframe with the similarity scores
    Taxa_stats = pd.DataFrame([(col, sum(df[col] > 0 ), similarity_scores[col]) for col in df.columns],
                              columns=['Col_name', 'Num_of_species', 'Similarity_jac']).set_index('Col_name')
    return Taxa_stats, ground_truth

def compute_bray_curtis_dissimilarity(df,Taxa_stats, ground_truth):
    '''Computes the Bray-Curtis dissimilarity between all columns in df and a specific column.'''
    def distance_calculation(u, v):
        l2 = ((u - v) ** 2).sum() ** 0.5
        l1 = sum(abs(u - v))
        den = sum(u + v)
        # distance = l1 / den
        return l1, l2

    # get the column to compare to
    u = df[ground_truth]
    # initialize an empty list to hold the dissimilarity scores
    dissimilarity_scores = []
    l1_scores = []
    l2_scores = []
    # loop over all columns
    for col in df.columns:
        # compute the Bray-Curtis dissimilarity between the two columns
        v = df[col]
        l1, l2 = distance_calculation(u, v)
        # add the distance score to the list of dissimilarity scores
        distance = braycurtis(u, v)
        dissimilarity_scores.append(distance)
        l1_scores.append(l1)
        l2_scores.append(l2)
    # create a DataFrame with the dissimilarity scores and return it
    dissimilarity_matrix = pd.DataFrame({'Dissimilarity': dissimilarity_scores, 'l1_score': l1_scores, 'l2_score': l2_scores},
                                        index=df.columns)
    df = pd.concat([Taxa_stats, dissimilarity_matrix], axis=1).reset_index()

    def get_sort_value(name):
        # return float(name.split('.')[1][:-1])
        return name

    stat_species = df.sort_values(by=df.columns[0], key=lambda x: x.apply(get_sort_value), ascending=True).set_index(df.columns[0]).round(2)
    return stat_species

def min_max_species(df,stat_species):
    stat_species['Min_Abundance'] = df.iloc[:, :-1].min(axis=0)
    stat_species['Max_Abundance'] = df.iloc[:, :-1].max(axis=0)
    stat_species.drop(stat_species.index[0], inplace=True)
    print(stat_species)

def main():

    if len(sys.argv) != 3:
        quit("\nUsage: " + sys.argv[0] + " <path_fastq_file_X_cov> <Taxa_list_ref> \n\n")

    path_fastq_file_X_cov = sys.argv[1]
    Taxa_list_ref = sys.argv[2]

    file_paths, id_files = Get_id(path_fastq_file_X_cov)
    Taxa_list_ref = Create_ref(Taxa_list_ref)

    if not os.path.isdir("mp4.Simulation"):
        os.mkdir("mp4.simulation")
    for path_file, id_file in zip(file_paths, id_files):
        metaphlan_file_out = f"{id_file}.txt"
        os.system(f'/data1/software/metaphlan/run-metaphlan.sh {path_file} ~/metaanalysis/mp4.simulation/{metaphlan_file_out}'
                  f' 40 > ~/metaanalysis/mp4.simulation/{metaphlan_file_out}.stdout')

    Taxa_sim = Mp4_txt_2_csv()
    combine2ref = cal_dis(Taxa_sim, Taxa_list_ref)
    print(Taxa_sim, Taxa_list_ref)
    Taxa_stats, ground_truth = compute_jaccard_scores(combine2ref)
    stat_species = compute_bray_curtis_dissimilarity(combine2ref,Taxa_stats, ground_truth)
    min_max_species(combine2ref, stat_species)

if __name__ == '__main__':
    main()