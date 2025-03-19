#!/usr/bin/python3
# Upload packages
from skbio.diversity import beta_diversity, alpha
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from fancyimpute import IterativeImputer
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from skbio.stats.composition import ancom, multiplicative_replacement
import itertools
from statsmodels.stats.power import TTestIndPower
from sklearn.metrics import roc_auc_score
from scipy.spatial import cKDTree
from tableone import TableOne
from IPython.display import display
import os
import warnings
from skbio import DistanceMatrix
from statannot import add_stat_annotation
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


# Suppress warnings
warnings.filterwarnings('ignore')


class MicrobiomeAnalysis:
    def __init__(self):
        # Initialize datasets to None
        self.IBD_1000 = None
        self.Eran_Elinav = None
        #self.Gevers = None
        self.HMP1 = None
        self.HMP2 = None
        self.HMP3 = None
        self.ibdmdb = None
        self.ibdmdb2 = None
        self.Lewis = None
        self.MetaHit = None
        self.MetaHit2 = None
        self.disease_counts = None
        self.LLDeep_PRIZM = None
        self.LLDeep = None
        self.PRISM = None 
        self.color_palette = {
            "1000 IBD": "#001219",
            "Federici et al., (2022)": "#005f73",
            "HMP1": "#0a9396",
            "HMP2": "#83c5be",
            "HMP2 Pilot": "#80B1D3",
            "HMP3": "#ee9b00",
            "IBDMDB": "#ca6702",
            "IBDMDB2": "#bb3e03",
            "Lewis (2015)": "#ae2012",
            "MetaHit": "#9b2226",
            "MetaHit2": "#e9d8a6",
            "LLDeep": "#335c67",
            "PRISM": "#38a3a5",
            'Israel':"#ae2012",
            'France':"#bb3e03",
            'US':"#80B1D3",
            'Germany':"#005f73",
            'A':"#005f73",
            'B': '#ae2012'
             }

        n = 21
        for i in range(1, n + 1):
            # Generating a hexadecimal component for the color code
            hex_component_A = format((50 + i) % 256, '02x')  # For group A
            hex_component_B = format((150 + i) % 256, '02x')  # For group B

            color_code_A = f"#00{hex_component_A}73"
            color_code_B = f"#00{hex_component_B}73"

            self.color_palette[f'A_{i}'] = color_code_A
            self.color_palette[f'B_{i}'] = color_code_B


    def find_string(self, df, your_string):
        string_columns = df.select_dtypes(include='object')
        result = string_columns.apply(lambda x: x.str.contains(your_string)).any()
        return result.index[result]

    def search_columns(self, df, substring):
        return [col for col in df.columns if substring in col]

    def search_index(self, df: pd.DataFrame, search_term: str) -> pd.DataFrame:
        matching_indices = [index for index in df.index if search_term in str(index)]
        return df.loc[matching_indices]

    def read_mp4(self, path_mp4_files):
        dataframes = {}
        samples = [x for x in os.listdir(path_mp4_files) if x.endswith(".csv")]
        for sample in samples:
            csv_name = sample.split('_species_abundance.csv')[0]
            df = pd.read_csv(f'{path_mp4_files}/{sample}')
            df = df.set_index(df.columns[0]).T
            df = df.loc[:, ~df.columns.str.startswith('t__')]
            dataframes[csv_name] = df
        sorted_keys = sorted(dataframes.keys())
        return dataframes, sorted_keys
    
    def pval_to_stars(self, pval):
        if pval > 0.05:
            return 'ns'
        elif pval > 0.01:
            return '*'
        elif pval > 0.001:
            return '**'
        elif pval > 0.0001:
            return '***'
        else:
            return '****'

    
    def compute_beta_diversities(self, df_list, metric, TAXA_level, threshold):
        dfs = []
        medians_dict = {}
        df_dict_projects = {}
        df_names = []

        for df_name, df in df_list:
            df = df.apply(pd.to_numeric, errors='coerce')
            df = df.where(df >= threshold)
            df = df.fillna(0)
            df = df.loc[:, (df != 0).any(axis=0)]  # drop columns if all abundance contains 0

            if TAXA_level == "g":
                df = df.loc[:, df.columns.str.startswith('g__')]
            elif TAXA_level == "s":
                df = df.loc[:, df.columns.str.startswith('s__')]

            ids = df.index.tolist()

            data = df.values
            # Convert to binary if metric is Jaccard
            if metric == "jaccard":
                data = np.where(data != 0, 1, 0) 
            dm = beta_diversity(metric, data, ids)
            df_dict_projects[df_name] = dm

            distances_df = pd.DataFrame(dm.condensed_form(), columns=[df_name])
            dfs.append(distances_df)

            median_value = distances_df[df_name].median()
            if metric == "jaccard":
                medians_dict[df_name] = 1-median_value
            else:
                medians_dict[df_name] = median_value
            df_names.append(df_name)

        all_pairs = list(itertools.combinations(df_names, 2))

        distances_df = pd.concat(dfs, axis=1)
        if metric == "jaccard":
            distances_df = 1 - distances_df

        sns.set_palette(list(self.color_palette.values()))
        plt.figure(figsize=(10, 8), dpi=800)
        ax = sns.boxplot(data=distances_df, palette=self.color_palette, linewidth=0.5)
        labels = ax.get_xticklabels()
        ax.set_xticklabels(labels, rotation=90)
        plt.tight_layout()  # Adjust layout to ensure labels fit

        valid_pairs_beta = [pair for pair in all_pairs if pair[0] in distances_df.columns and pair[1] in distances_df.columns]
        ax,test_results =  add_stat_annotation(ax, data=distances_df, box_pairs=valid_pairs_beta, test='Mann-Whitney', text_format='star', loc='outside', verbose=2)

        if TAXA_level == 'g':
            title_taxa = 'Genus'
        elif TAXA_level == 's':
            title_taxa = 'Species'
        else:
            title_taxa = "OTU"  

        # Function to convert p-values to significance stars
        def pval_to_stars(pval):
            if pval > 0.05:
                return ' '
            elif pval > 0.01:
                return '*'
            elif pval > 0.001:
                return '**'
            elif pval > 0.0001:
                return '***'
            else:
                return '***'

        # Extract dataset names from the list of valid pairs
        datasets = set()
        for pair in valid_pairs_beta:
            datasets.update(pair)

        # Initialize an empty DataFrame with the datasets as indices and columns
        pval_matrix = pd.DataFrame(index=sorted(datasets), columns=sorted(datasets))

        # Fill in the diagonal with empty value since it represents the same dataset comparison
        np.fill_diagonal(pval_matrix.values, ' ')

        # Iterate through each StatResult object in test_results
        for test_result in test_results:
            # Extract the dataset names from the comparison
            dataset1, dataset2 = test_result.box1, test_result.box2
            # Fill in the corresponding cell with the p-value or statistical significance
            pval_matrix.loc[dataset1, dataset2] = test_result.pval
            pval_matrix.loc[dataset2, dataset1] = test_result.pval  # Ensure symmetry

        # Replace p-values with significance stars
        self.pval_matrix_stars = pval_matrix.applymap(lambda x: pval_to_stars(x) if isinstance(x, float) else x)

        plt.title(f'Median {metric.capitalize()} β-diversity Between Each Pair of {title_taxa} Tables')
        plt.show()
        return medians_dict, df_dict_projects


    def compute_median_diversities(self, df_list, metric, TAXA_level,threshold):
        median_diversities = pd.DataFrame(np.zeros((len(df_list), len(df_list))), columns=[x[0] for x in df_list],
                                          index=[x[0] for x in df_list])

        for i in range(len(df_list)):
            for j in range(i + 1, len(df_list)):

                df_name_i, df_i = df_list[i]
                df_name_j, df_j = df_list[j]

                df_i = df_i.apply(pd.to_numeric, errors='coerce')
                df_i = df_i.where(df_i >= threshold)
                df_i = df_i.fillna(0)
                df_i = df_i.loc[:, (df_i != 0).any(axis=0)]
                df_i = df_i.sort_index(axis=1).astype(int)
                if TAXA_level == "g":
                    df_i = df_i.loc[:, df_i.columns.str.startswith('g__')]
                    data_i = df_i.values
                elif TAXA_level == "s":
                    df_i = df_i.loc[:, df_i.columns.str.startswith('s__')]
                    data_i = df_i.values

                df_j = df_j.apply(pd.to_numeric, errors='coerce')
                df_j = df_j.where(df_j >= threshold)
                df_j = df_j.fillna(0)
                df_j = df_j.loc[:, (df_j != 0).any(axis=0)]
                df_j = df_j.sort_index(axis=1).astype(int)
                if TAXA_level == "g":
                    df_j = df_j.loc[:, df_j.columns.str.startswith('g__')]
                    data_j = df_j.values
                elif TAXA_level == "s":
                    df_j = df_j.loc[:, df_j.columns.str.startswith('s__')]
                    data_j = df_j.values
                    
                if metric == "jaccard":
                    df_i = df_i.where(df_i == 0, 1)
                    df_j = df_j.where(df_j == 0, 1)

                # Align the dataframes
                df_i, df_j = df_i.align(df_j, axis=1, fill_value=0)

                # Now, convert them to NumPy arrays for diversity calculation
                data_i = df_i.values
                data_j = df_j.values
                
                # Compute Bray-Curtis beta diversity
                dm_between = beta_diversity(metric, np.concatenate([data_i, data_j]),
                                            ids=[f"{df_name_i}_{id_}" for id_ in df_i.index]
                                                + [f"{df_name_j}_{id_}" for id_ in df_j.index])

                dm = pd.DataFrame(dm_between.data, index=dm_between.ids, columns=dm_between.ids)
                dm = dm.fillna(0)
                median_diversity = np.median(dm)
                median_diversities.loc[df_name_i, df_name_j] = median_diversity
                median_diversities.loc[df_name_j, df_name_i] = median_diversity
        if metric == "jaccard":
            median_diversities = 1 - median_diversities
                    
        medians_dict,_ = self.compute_beta_diversities(df_list, metric, TAXA_level,threshold)
        for index, (key, value) in enumerate(medians_dict.items()):
            if key in median_diversities.index and key in median_diversities.columns:
                median_diversities.at[key, key] = value
        
        plt.figure(figsize=(10, 8), dpi=800)

        if TAXA_level == 'g':
            taxa_label = 'genus'
        elif TAXA_level == 's':
            taxa_label = 'species'
        else:
            taxa_label = 'OTU'
            
        median_diversities_rounded = median_diversities.round(2)
        # Assuming self.pval_matrix_stars is a DataFrame with the same indices and columns as median_diversities
        # and contains strings with significance stars like "*", "**", "***", etc.
        display(median_diversities_rounded)
        # Add a newline character before each significance star in the annotations
        median_diversities_rounded = median_diversities_rounded.reindex_like(self.pval_matrix_stars)
        self.pval_matrix_stars = self.pval_matrix_stars.reindex_like(median_diversities_rounded)
        annotations_with_newline = median_diversities_rounded.astype(str) #+ '\n' + self.pval_matrix_stars.fillna('')
        #sns.set_palette(list(self.color_palette.values()))
        display(annotations_with_newline)
        # Use the numeric DataFrame for the heatmap data and the annotations DataFrame for annotating
        sns.clustermap(median_diversities_rounded, annot=annotations_with_newline, fmt='', cmap=sns.color_palette("icefire", as_cmap=True))

        plt.title(f'Median {metric.capitalize()} β-diversity Between Each Pair of {taxa_label} Tables')
        plt.show()
        return medians_dict
        
    def compute_alpha_diversities(self, df_list, TAXA_level, threshold):
        shannon_diversities = []
        chao1_diversities = []
        simpson_diversities = []
        df_names = []

        # Loop through the DataFrame list
        for df_name, df in df_list:
            df = df.apply(pd.to_numeric, errors='coerce')
            df = df.where(df >= threshold)
            df = df.fillna(0)
            df = df.loc[:, (df != 0).any(axis=0)]  # drop columns if all abundance contains 0

            if TAXA_level == "g":
                df = df.loc[:, df.columns.str.startswith('g__')]
            elif TAXA_level == "s":
                df = df.loc[:, df.columns.str.startswith('s__')]

            # Check if dataframe is not empty after filtering
            if not df.empty:
                df_names.append(df_name)
                # Compute alpha diversities
                shannon_diversities.append(df.apply(alpha.shannon, axis=1))
                chao1_diversities.append(df.apply(alpha.chao1, axis=1))
                simpson_diversities.append(df.apply(alpha.simpson, axis=1))

        # Prepare the results for plotting
        results_df_shannon = pd.concat(shannon_diversities, axis=1, keys=df_names)
        results_df_chao1 = pd.concat(chao1_diversities, axis=1, keys=df_names)
        results_df_simpson = pd.concat(simpson_diversities, axis=1, keys=df_names)

        # Define the list of projects
        all_pairs = list(itertools.combinations(df_names, 2))

        sns.set_palette(list(self.color_palette.values()))

        # Plot results
        plt.figure(figsize=(15, 8), dpi=800)

        ax1 = plt.subplot(1, 3, 1)
        sns.boxplot(data=results_df_shannon, palette=self.color_palette, linewidth=0.5)
        plt.title('Shannon Diversity')
        plt.xticks(rotation=90)
        valid_pairs_shannon = [pair for pair in all_pairs if pair[0] in results_df_shannon.columns and pair[1] in results_df_shannon.columns]
        add_stat_annotation(ax1, data=results_df_shannon, box_pairs=valid_pairs_shannon, test='Mann-Whitney', text_format='star', loc='outside', verbose=2)

        # Chao1
        ax2 = plt.subplot(1, 3, 2)
        sns.boxplot(data=results_df_chao1, palette=self.color_palette, linewidth=0.5)
        plt.title('Chao1 Diversity')
        plt.xticks(rotation=90)
        valid_pairs_chao1 = [pair for pair in all_pairs if pair[0] in results_df_chao1.columns and pair[1] in results_df_chao1.columns]
        add_stat_annotation(ax2, data=results_df_chao1, box_pairs=valid_pairs_chao1, test='Mann-Whitney', text_format='star', loc='outside', verbose=2)

        # Simpson
        ax3 = plt.subplot(1, 3, 3)
        sns.boxplot(data=results_df_simpson, palette=self.color_palette, linewidth=0.5)
        plt.title('Simpson Diversity')
        plt.xticks(rotation=90)
        valid_pairs_simpson = [pair for pair in all_pairs if pair[0] in results_df_simpson.columns and pair[1] in results_df_simpson.columns]
        add_stat_annotation(ax3, data=results_df_simpson, box_pairs=valid_pairs_simpson, test='Mann-Whitney', text_format='star', loc='outside', verbose=2)

        plt.tight_layout()
        plt.show()

    def propagate_values(self, df, groupby_col, target_cols):
        for target_col in target_cols:
            missing_values = df[df[target_col].isna()][groupby_col]
            existing_values = df.dropna(subset=[target_col]).groupby(groupby_col)[target_col].first()
            df.loc[df[groupby_col].isin(missing_values), target_col] = df[groupby_col].map(existing_values)
        return df

    def process_data(self):
        # Load datasets
        self.IBD_1000 = pd.read_csv(
            "/data1/Human/1000_IBD_cohort/1000_IBD_cohort.clinincal_data.EGAD00001003991/EGA_metagenomics_Phenotypes_1000IBD_release_2.txt",
            delimiter="\t").set_index("ID_1000IBD")
        self.Eran_Elinav = pd.read_csv("/home/odeds/metaanalysis/datasets/multi-cohort_dataset/Eran_Elinav/metadata.txt",
                                       delimiter="\t").set_index("biosample")

        self.HMP1 = pd.read_csv("/home/odeds/metaanalysis/datasets/multi-cohort_dataset/HMP1/Stool/HMP.phase1.info.txt", delimiter="\t")
        self.HMP2 = pd.read_csv("/home/odeds/metaanalysis/datasets/multi-cohort_dataset/HMP2/Stool/HMP.phase2.info.txt", delimiter="\t",
                           skiprows=1)
        self.HMP3 = pd.read_csv("/home/odeds/metaanalysis/datasets/multi-cohort_dataset/HMP3/Stool/HMP.phase3.info.txt", delimiter="\t")
        self.ibdmdb = pd.read_csv("/home/odeds/metaanalysis/datasets/multi-cohort_dataset/ibdmdb/info.txt", delimiter="\t", skiprows=1).set_index(
            "External ID")
        self.ibdmdb = self.ibdmdb[self.ibdmdb["data_type"] == "metagenomics"]
        self.ibdmdb = analysis.propagate_values(self.ibdmdb, 'Participant ID', ["specify_race"])

        self.ibdmdb2 = pd.read_csv("/home/odeds/metaanalysis/datasets/multi-cohort_dataset/ibdmdb2/hmp2_metadata.csv")
        self.ibdmdb2 = analysis.propagate_values(self.ibdmdb2, 'Participant ID', ['BMI', "Specify race"])

        self.HMP2_pilot = self.ibdmdb2[(self.ibdmdb2["data_type"] == "metagenomics") & (self.ibdmdb2["External ID"].str.endswith("_P"))]
        self.ibdmdb2 = self.ibdmdb2[(self.ibdmdb2["data_type"] == "metagenomics") & (~self.ibdmdb2["External ID"].str.endswith("_P"))]

        self.ibdmdb2.set_index("External ID", inplace=True)
        self.HMP2_pilot.set_index("External ID", inplace=True)

        self.Lewis = pd.read_csv("/home/odeds/metaanalysis/datasets/multi-cohort_dataset/Lewis/SRP057027.info.txt", delimiter="\t").set_index(
            "sample_accession")
        self.MetaHit = pd.read_csv("/home/odeds/metaanalysis/datasets/multi-cohort_dataset/MetaHit/Study.design.tab", delimiter="\t").set_index(
            "Sample")
        self.MetaHit2 = pd.read_csv("/home/odeds/metaanalysis/datasets/multi-cohort_dataset/MetaHit2/nbt.2939-S2.with-sample-IDs.txt",
                               delimiter="\t").set_index(
            "Sample")

        self.LLDeep_PRIZM = pd.read_csv("~/metaanalysis/datasets/multi-cohort_dataset/LLDeep_PRIZM/metadata2.tsv",
                               delimiter="\t").set_index("Unnamed: 0")
        self.LLDeep = self.LLDeep_PRIZM[self.LLDeep_PRIZM["Sample"].str.startswith("LLDeep")]
        self.PRISM = self.LLDeep_PRIZM[self.LLDeep_PRIZM["Sample"].str.startswith("PRISM")]


        #This section shows the number of patients by cohort distribution—time series or not.
        cohorts_non_ts = [self.HMP1, self.HMP2, self.HMP3, self.MetaHit, self.IBD_1000, self.Eran_Elinav, self.LLDeep,
                          self.PRISM]
        cohorts_ts = [self.ibdmdb, self.ibdmdb2, self.HMP2_pilot, self.Lewis, self.MetaHit2]

        sum_non_ts = sum(df.shape[0] for df in cohorts_non_ts)
        sum_ts = sum(df.shape[0] for df in cohorts_ts)

        print(f'sum of patients from cohorts (not TS): {sum_non_ts}')
        print(f'sum of patients from Time series (TS) cohorts: {sum_ts}')
        print("#########################################################")
        
        
        # Process the data
        # this output from Sort_Metaphlan.py (into same directory)
        dfs, sorted_keys = self.read_mp4("/home/odeds/metaanalysis/datasets/species_abundance/")

        dfs['1000_IBD_UC']['Diagnosis'] = 'UC'
        dfs['1000_IBD_CD']['Diagnosis'] = 'CD'

        dfs['ibdmdb_nonIBD']['Diagnosis'] = 'nonIBD'
        dfs['ibdmdb_UC']['Diagnosis'] = 'UC'
        dfs['ibdmdb_CD']['Diagnosis'] = 'CD'

        dfs['ibdmdb2_nonIBD']['Diagnosis'] = 'nonIBD'
        dfs['ibdmdb2_CD']['Diagnosis'] = 'CD'
        dfs['ibdmdb2_UC']['Diagnosis'] = 'UC'

        self.IBD_1000_T = pd.concat([dfs['1000_IBD_UC'], dfs['1000_IBD_CD'], dfs['1000_IBD_IBDU']], axis=0)

        self.ibdmdb_T = pd.concat([dfs['ibdmdb_CD'], dfs['ibdmdb_UC'], dfs['ibdmdb_nonIBD']], axis=0)

        self.ibdmdb2_T = pd.concat([dfs['ibdmdb2_CD'], dfs['ibdmdb2_UC'], dfs['ibdmdb2_nonIBD']], axis=0)

        self.Eran_Elinav_T = dfs['Eran_Elinav']

        self.HMP1_T = dfs['HMP1']

        self.HMP2_T = dfs['HMP2']

        self.HMP2_pilot_T = dfs['HMP2_pilot']

        self.HMP3_T = dfs['HMP3']

        self.Lewis_T = dfs['Lewis']

        self.MetaHit_T = dfs['MetaHit']
        
        self.MetaHit2_T = dfs['MetaHit2']

        self.LLDeep_T = dfs['LLDeep']

        self.PRISM_T = dfs['PRISM']
        
         #This section shows the number of patients by cohort distribution—time series or not.
        cohorts_non_ts = [self.HMP1_T, self.HMP2_T, self.HMP3_T, self.MetaHit_T, self.IBD_1000_T, self.Eran_Elinav_T, self.LLDeep_T,self.PRISM_T]
        cohorts_ts = [self.ibdmdb_T, self.ibdmdb2_T, self.HMP2_pilot_T, self.Lewis_T, self.MetaHit2_T]

        sum_non_ts = sum(df.shape[0] for df in cohorts_non_ts)
        sum_ts = sum(df.shape[0] for df in cohorts_ts)

        print(f'sum of samples from cohorts (not TS): {sum_non_ts}')
        print(f'sum of samples from Time series (TS) cohorts: {sum_ts}')
        print("#########################################################")

        indices_shared = self.IBD_1000.index.intersection(self.IBD_1000_T.index)
        self.IBD_1000 = self.IBD_1000.loc[indices_shared].reset_index()

        indices_shared = self.ibdmdb.index.intersection(self.ibdmdb_T.index)
        self.ibdmdb = self.ibdmdb.loc[indices_shared].reset_index()

        indices_shared = self.ibdmdb2.index.intersection(self.ibdmdb2_T.index)
        self.ibdmdb2 = self.ibdmdb2.loc[indices_shared].reset_index()

        indices_shared = self.Eran_Elinav.index.intersection(self.Eran_Elinav_T.index)
        self.Eran_Elinav = self.Eran_Elinav.loc[indices_shared].reset_index()

        indices_shared = self.HMP1.index.intersection(self.HMP1_T.index)
        self.HMP1 = self.HMP1.loc[indices_shared].reset_index()

        indices_shared = self.HMP2.index.intersection(self.HMP2_T.index)
        self.HMP2 = self.HMP2.loc[indices_shared].reset_index()

        indices_shared = self.HMP2_pilot.index.intersection(self.HMP2_pilot_T.index)
        self.HMP2_pilot = self.HMP2_pilot.loc[indices_shared].reset_index()

        indices_shared = self.HMP3.index.intersection(self.HMP3_T.index)
        self.HMP3 = self.HMP3.loc[indices_shared].reset_index()

        indices_shared = self.Lewis.index.intersection(self.Lewis_T.index)
        self.Lewis = self.Lewis.loc[indices_shared].reset_index()

        indices_shared = self.MetaHit.index.intersection(self.MetaHit_T.index)
        self.MetaHit = self.MetaHit.loc[indices_shared].reset_index()
        
        indices_shared = self.MetaHit2.index.intersection(self.MetaHit2_T.index)
        self.MetaHit2 = self.MetaHit2.loc[indices_shared].reset_index()

        indices_shared = self.LLDeep.index.intersection(self.LLDeep_T.index)
        self.LLDeep = self.LLDeep.loc[indices_shared].reset_index()

        indices_shared = self.PRISM.index.intersection(self.PRISM_T.index)
        self.PRISM = self.PRISM.loc[indices_shared].reset_index()

        # Befor phase 0 
        cohorts_non_ts = [self.HMP1, self.HMP2, self.HMP3, self.MetaHit, self.IBD_1000, self.Eran_Elinav, self.LLDeep,
                          self.PRISM]
        cohorts_ts = [self.ibdmdb, self.ibdmdb2, self.HMP2_pilot, self.Lewis, self.MetaHit2]

        sum_non_ts = sum(df.shape[0] for df in cohorts_non_ts)
        sum_ts = sum(df.shape[0] for df in cohorts_ts)

        print(f'Befor phase 0  -  sum of patients from cohorts (not TS): {sum_non_ts}')
        print(f'Befor phase 0  - sum of patients from Time series (TS) cohorts: {sum_ts}')
        print("#########################################################")
        
        
         #This section shows the number of patients by cohort distribution—time series or not.
        cohorts_non_ts = [self.HMP1_T, self.HMP2_T, self.HMP3_T, self.MetaHit_T, self.IBD_1000_T, self.Eran_Elinav_T, self.LLDeep_T,self.PRISM_T]
        cohorts_ts = [self.ibdmdb_T, self.ibdmdb2_T, self.HMP2_pilot_T, self.Lewis_T, self.MetaHit2_T]

        sum_non_ts = sum(df.shape[0] for df in cohorts_non_ts)
        sum_ts = sum(df.shape[0] for df in cohorts_ts)

        print(f'sum of samples from cohorts (not TS): {sum_non_ts}')
        print(f'sum of samples from Time series (TS) cohorts: {sum_ts}')
        print("#########################################################")
        # - Time series data structures are dealt with in the following project cell.
        # - One row is generated at random for each patient.
        self.ibdmdb = self.ibdmdb.groupby('Participant ID').first()
        self.HMP2_pilot = self.HMP2_pilot.groupby('Participant ID').first()
        self.ibdmdb2 = self.ibdmdb2.groupby('Participant ID').first()
        self.Lewis['diagnosis'] = 'CD'
        self.Lewis.loc[:25, 'diagnosis'] = np.where(self.Lewis.iloc[:26].isna().any(axis=1), 'nonIBD', 'CD')
        self.Lewis = self.Lewis.groupby('Subject').first()
        self.MetaHit2 = self.MetaHit2.groupby('Individual ID').first()
        self.MetaHit = self.MetaHit.groupby('Name').first()

        # Summarizing diseases
        index = ['CD', 'UC', 'IBD', 'nonIBD']
        self.disease_counts = pd.DataFrame(index=index)

        self.IBD_1000['diagnosis_last_record'] = self.IBD_1000['diagnosis_last_record'].replace({"IBDU": "IBD"})
        self.disease_counts['IBD_1000'] = self.IBD_1000['diagnosis_last_record'].value_counts()

        self.Eran_Elinav['Disease'] = self.Eran_Elinav['Disease'].replace({"Ctrl": "nonIBD"})
        self.disease_counts['EranElinav'] = self.Eran_Elinav['Disease'].value_counts()

#         self.['Diagnosis'] = self.Gevers['Diagnosis'].replace({"Not IBD": "nonIBD", "IBD-U": "IBD"})
#         self.disease_counts['Gevers'] = self.Gevers['Diagnosis'].value_counts()

        self.HMP1["study design"] = self.HMP1["study design"].replace({' Control Set': "nonIBD"})
        self.disease_counts['HMP1'] = self.HMP1["study design"].value_counts()

        self.HMP2["study design"] = self.HMP2["study design"].replace({' Control Set': "nonIBD"})
        self.disease_counts['HMP2'] = self.HMP2["study design"].value_counts()

        indicate = (self.HMP3["isolation_source"] == 'GUT').sum()
        self.disease_counts.at['nonIBD', 'HMP3'] = indicate

        self.ibdmdb["diagnosis"] = self.ibdmdb["diagnosis"].replace(
            {"Crohn's Disease": "CD", "Ulcerative colitis": "UC", "Healthy control": "nonIBD"})
        self.disease_counts["IBDMDB"] = self.ibdmdb["diagnosis"].value_counts()

        self.disease_counts["HMP2_pilot"] = self.HMP2_pilot["diagnosis"].value_counts()

        self.disease_counts["IBDMDB2"] = self.ibdmdb2["diagnosis"].value_counts()

        self.disease_counts["Lewis"] = self.Lewis['diagnosis'].value_counts()

        condition1 = self.MetaHit['IBD'] == 'N'
        condition2 = (self.MetaHit['IBD'] == 'Y') & (self.MetaHit.index.str.contains('UC'))
        condition3 = (self.MetaHit['IBD'] == 'Y') & (self.MetaHit.index.str.contains('CD'))
        self.MetaHit['diagnosis'] = np.where(condition1, 'nonIBD', np.where(condition2, 'UC', np.where(condition3, 'CD', 'IBD')))

        self.disease_counts["MetaHit"] = self.MetaHit['diagnosis'].value_counts()
        
        self.MetaHit2["Health Status"] = self.MetaHit2["Health Status"].replace(
            {"Crohns disease": "CD", "Ulcerative colitis": "UC", "Healthy": "nonIBD"})    
        self.disease_counts['MetaHit2'] = self.MetaHit2["Health Status"].value_counts()

        self.LLDeep["Study.Group"] = self.LLDeep["Study.Group"].replace({"Control": "nonIBD"})
        self.disease_counts['LLDeep'] = self.LLDeep["Study.Group"].value_counts()

        self.PRISM["Study.Group"] = self.PRISM["Study.Group"].replace({"Control": "nonIBD"})
        self.disease_counts['PRISM'] = self.PRISM["Study.Group"].value_counts()

        self.disease_counts = self.disease_counts.fillna(0).astype(int)
        self.disease_counts['Total_Patients'] = self.disease_counts.apply(sum, axis=1)
        self.disease_counts.index.name = 'Diagnosis'
        
        
        #This section shows the number of patients by cohort distribution—time series or not.
        cohorts_non_ts = [self.HMP1, self.HMP2, self.HMP3, self.MetaHit, self.IBD_1000, self.Eran_Elinav, self.LLDeep,
                          self.PRISM]
        cohorts_ts = [self.ibdmdb, self.ibdmdb2, self.HMP2_pilot, self.Lewis, self.MetaHit2]

        sum_non_ts = sum(df.shape[0] for df in cohorts_non_ts)
        sum_ts = sum(df.shape[0] for df in cohorts_ts)

        print(f'After cross -  sum of patients from cohorts (not TS): {sum_non_ts}')
        print(f'After cross - sum of patients from Time series (TS) cohorts: {sum_ts}')
        print("#########################################################")
         #This section shows the number of patients by cohort distribution—time series or not.
        cohorts_non_ts = [self.HMP1_T, self.HMP2_T, self.HMP3_T, self.MetaHit_T, self.IBD_1000_T, self.Eran_Elinav_T, self.LLDeep_T,self.PRISM_T]
        cohorts_ts = [self.ibdmdb_T, self.ibdmdb2_T, self.HMP2_pilot_T, self.Lewis_T, self.MetaHit2_T]

        sum_non_ts = sum(df.shape[0] for df in cohorts_non_ts)
        sum_ts = sum(df.shape[0] for df in cohorts_ts)

        print(f'sum of samples from cohorts (not TS): {sum_non_ts}')
        print(f'sum of samples from Time series (TS) cohorts: {sum_ts}')
        print("#########################################################")
        return {
            'IBD_1000': self.IBD_1000,
            'IBD_1000_T': self.IBD_1000_T,
            'Eran_Elinav': self.Eran_Elinav,
            'Eran_Elinav_T': self.Eran_Elinav_T,
            #'Gevers': self.Gevers,
            #'Gevers_T': self.Gevers_T,
            'HMP1': self.HMP1,
            'HMP1_T': self.HMP1_T,
            'HMP2': self.HMP2,
            'HMP2_T': self.HMP2_T,
            'HMP2_pilot': self.HMP2_pilot,
            'HMP2_pilot_T': self.HMP2_pilot_T,
            'HMP3': self.HMP3,
            'HMP3_T': self.HMP3_T,
            'ibdmdb': self.ibdmdb,
            'ibdmdb_T': self.ibdmdb_T,
            'ibdmdb2': self.ibdmdb2,
            'ibdmdb2_T': self.ibdmdb2_T,
            'Lewis': self.Lewis,
            'Lewis_T': self.Lewis_T,
            'MetaHit': self.MetaHit,
            'MetaHit_T': self.MetaHit_T,
            'MetaHit2': self.MetaHit2,
            'MetaHit2_T': self.MetaHit2_T,
            'LLDeep': self.LLDeep,
            'LLDeep_T': self.LLDeep_T,
            'PRISM': self.PRISM,
            'PRISM_T': self.PRISM_T,
            'disease_counts': self.disease_counts,
            'sorted_keys': sorted_keys
        }
# To use the class:
analysis = MicrobiomeAnalysis()
results = analysis.process_data()
