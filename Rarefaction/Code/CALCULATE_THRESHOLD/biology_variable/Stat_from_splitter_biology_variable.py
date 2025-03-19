
#!/usr/bin/python3

# Import packages
import os
import pandas as pd
import sys

class MetaphlanAnalysis_from_shuffle:
    def __init__(self, path_mp4_files, directory_name):
        self.path_mp4_files = path_mp4_files
        self.directory_name = directory_name

    def get_data(self, threshold):
        species_main = None
        samples = [os.path.join(dp, f) for dp, dn, filenames in os.walk(self.path_mp4_files) for f in filenames if
                   f.endswith('.txt')]
        for sample in samples:
            species_abundance = dict()
            self.n_reads = sample.split(".")[1].split("_")[1].split(".txt")[0]
            sample_name = sample.split("/")[-1].split(".txt")[0]
            with open(sample) as s:
                f = s.readlines()
                for line in f:
                    if line.startswith("k__Bacteria"):
                        try:
                            species = line.strip().split("\t")[0].split("|")[-1]
                            abundance = float(line.split()[2])
                            if species.startswith('s__') and abundance >= threshold:
                                species_abundance[species] = abundance
                        except ValueError:
                            print(f'This value: {line.split()[2]} is not an abundance')
            Taxa_Abundance = pd.DataFrame(species_abundance.items(), columns=['Species', sample_name]).set_index('Species')
            species_main = pd.concat([species_main, Taxa_Abundance], axis=1, join="outer").fillna(0)
        return species_main


    def compute_jaccard_scores(self, df):
        similarity_scores = []
        cols = df.columns
        for i in range(len(cols)):
            for j in range(i + 1, len(cols)):
                col, feat = cols[i], cols[j]
                species1 = df.loc[df[col] > 0].index
                species2 = df.loc[df[feat] > 0].index
                m_species = (len(species1) + len(species2)) / 2
                intersection = len(set(species1).intersection(species2))
                union = len(set(species1).union(species2))
                score = intersection / union if union != 0 else 0
                feature_versusf = ' VS '.join(sorted([col, feat]))
                similarity_scores.append([feature_versusf, score, m_species])
        taxa_stats = pd.DataFrame(similarity_scores, columns=["index", "Jaccard Score", "m_species"]).set_index("index")
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

    def sort_df(self, df, threshold):
        df.reset_index(inplace=True)
        df.set_index("index", inplace=True)
        Taxa_stats = df.round(decimals=3)
        path = f'{self.directory_name}/Output_bio_MetaHit_versus_{threshold}.csv'
        Taxa_stats.to_csv(path)

    def run_analysis_from_shuffle(self):
        thresholds = [0, 0.001, 0.01, 0.1, 1]
        for threshold in thresholds:
            species_main = self.get_data(threshold)
            Taxa_stats = self.compute_jaccard_scores(species_main)
            df = self.compute_bray_curtis_dissimilarity(species_main, Taxa_stats)
            self.sort_df(df, threshold)

def main():
    if len(sys.argv) != 2:
        print(f"\nUsage: {sys.argv[0]} <Path to the Metaphlan output files>\n")
        sys.exit(1)

    path_mp4_files = sys.argv[1]

    # Assuming the output directory is the same as the input directory for Metaphlan files
    directory_name_out = os.path.dirname(path_mp4_files)

    analysis = MetaphlanAnalysis_from_shuffle(path_mp4_files, directory_name_out)
    analysis.run_analysis_from_shuffle()

if __name__ == '__main__':
    main()