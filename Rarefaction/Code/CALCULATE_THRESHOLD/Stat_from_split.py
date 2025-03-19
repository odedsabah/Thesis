#!/usr/bin/python3



'''Add number of species'''
# Upload packages
import os
import pandas as pd

'''
The script expects one command-line argument as input:

1. path_file_metaphlan_file: The path to the directory containing Metaphlan output files.

When executing the script from the command line, you would provide this argument. Here's an example:

<python stat_from_split.py /path/to/metaphlan_output_directory>

As for the output, the script generates a CSV file named "Taxa_stats.csv".
The file contains the sorted DataFrame with the calculated similarity and dissimilarity scores between pairs of samples.

The CSV file includes the following columns:
- Cols_name: The column names representing the pair of samples compared.
- jaccard_scores: The Jaccard similarity scores between the samples.
- Dissimilarity: The Bray-Curtis dissimilarity scores between the samples.
- l1_score: The L1 scores (sum of absolute differences) between the samples.
- l2_score: The L2 scores (Euclidean distances) between the samples.

The sorted DataFrame is saved to the CSV file in the same directory as the script.

Please note that the script assumes the presence of Metaphlan output files in the specified directory
,and it expects the files to be in a specific format. Ensure that the directory contains the necessary files for the script to execute successfully.
'''
class MetaphlanAnalysis_from_split:
    def __init__(self, path_mp4_files, directory_name):
        self.path_mp4_files = path_mp4_files
        self.directory_name = directory_name

    def get_data(self, threshold):
        species_main = None
        samples = [x for x in os.listdir(self.path_mp4_files) if x.endswith(".txt")]
        for sample in samples:
            species_abundance = dict()
            sample_name = sample[:-4]
            with open(f'{self.path_mp4_files}/{sample_name}.txt') as s:
                f = s.readlines()
                for line in f:
                    if line.startswith("k__Bacteria"):
                        try:
                            species = line.strip().split("\t")[0].split("|")[-1]
                            abundance = float(line.split()[2])
                            if species.startswith('s__') and abundance >= threshold:
                                species_abundance[species] = abundance
                        except:
                            print(f'this value: {line.split()[2]} is not an abundance')
            Taxa_Abundance = pd.DataFrame((species_abundance.items()), columns=['Species',sample_name]).set_index('Species')
            species_main = pd.concat([species_main, Taxa_Abundance], axis=1, join="outer").fillna(0)
        print(species_main)

        return species_main

    def compute_jaccard_scores(self, df):
        similarity_scores = {}
        cols = df.columns
        split_cols = [col.split("-") for col in cols]
        processed_cols = []

        for i in range(len(split_cols)):
            if cols[i] in processed_cols:
                continue
            for j in range(i + 1, len(split_cols)):
                if (split_cols[i][1].split('.')[0] == split_cols[j][1].split('.')[0] and split_cols[i][0] !=
                        split_cols[j][0]):
                    col, feat = cols[i], cols[j]
                    species1 = df.loc[df[col] > 0, :].index
                    species2 = df.loc[df[feat] > 0, :].index
                    m_species = (len(set(species1)) + len(set(species2))) / 2
                    intersection = len(set(species1).intersection(set(species2)))
                    union = len(set(species1).union(set(species2)))
                    score = intersection / union if union != 0 else 0
                    feature_versusf = ' VS '.join(sorted([col, feat]))
                    similarity_scores[feature_versusf] = [score, m_species]
                    processed_cols.append(cols[i])
                    processed_cols.append(cols[j])
                    break

        taxa_stats = pd.DataFrame(list(similarity_scores.items()), columns=['Cols_name', 'Score_and_m_species'])
        taxa_stats[['jaccard_scores', 'm_species']] = pd.DataFrame(taxa_stats['Score_and_m_species'].tolist(),
                                                                   index=taxa_stats.index)
        taxa_stats = taxa_stats.drop('Score_and_m_species', axis=1).set_index('Cols_name')

        return taxa_stats

    def compute_bray_curtis_dissimilarity(self, df, Taxa_stats):
        def distance_calculation(u, v):
            l2 = ((u - v) ** 2).sum() ** 0.5
            l1 = sum(abs(u - v))
            den = sum(u + v)
            distance = l1 / den if den != 0 else 0
            return distance, l1, l2

        dissimilarity_scores = {}
        l1_scores = {}
        l2_scores = {}

        for col_vs_feat in Taxa_stats.index:
            col, feat = col_vs_feat.split(' VS ')
            u = df[col]
            v = df[feat]
            distance, l1, l2 = distance_calculation(u, v)

            dissimilarity_scores[col_vs_feat] = distance
            l1_scores[col_vs_feat] = l1
            l2_scores[col_vs_feat] = l2

        dissimilarity_matrix = pd.DataFrame({'Dissimilarity': dissimilarity_scores, 'l1_score': l1_scores, 'l2_score': l2_scores})
        df = pd.concat([Taxa_stats, dissimilarity_matrix], axis=1)
        return df

    def parse_sort_key(self, item):
        parts = item.split(" VS ")
        keys = [(p.split('_')[-1].split('-')[0], p.split('_')[-1].split('-')[1].split('M')[0]) for p in parts]
        return min(k[1] for k in keys), min(k[0] for k in keys)

    def sort_df(self, df, threshold):
        df.reset_index(inplace=True)
        sorted_indices = sorted(df["index"], key=self.parse_sort_key)

        df.set_index("index", inplace=True)
        df = df.reindex(sorted_indices)
        Taxa_stats = df.round(decimals=3)
        id_sample = self.directory_name.split("/")[-1]
        path = f'{self.directory_name}/{id_sample}_split_{threshold}.csv'
        Taxa_stats.to_csv(path)

    def run_analysis_from_split(self, thresholds):
        for threshold in thresholds:
            species_main = self.get_data(threshold)
            Taxa_stats = self.compute_jaccard_scores(species_main)
            df = self.compute_bray_curtis_dissimilarity(species_main, Taxa_stats)
            self.sort_df(df, threshold)

