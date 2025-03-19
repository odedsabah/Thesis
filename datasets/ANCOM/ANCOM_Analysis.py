#!/usr/bin/python3

# Upload packages
from skbio.stats.composition import ancom, multiplicative_replacement
from Summary_Project import analysis, results
import itertools
import sys
import os
import pandas as pd

IBD_1000 = results["IBD_1000"]
IBD_1000_T = results["IBD_1000_T"]

Eran_Elinav = results['Eran_Elinav']
Eran_Elinav_T = results['Eran_Elinav_T']

Gevers = results['Gevers']
#Gevers_T = results['Gevers_T']

HMP1 = results['HMP1']
HMP1_T = results['HMP1_T']

HMP2 = results['HMP2']
HMP2_T = results['HMP2_T']

HMP2_pilot = results['HMP2_pilot']
HMP2_pilot_T = results['HMP2_pilot_T']

HMP3 = results['HMP3']
HMP3_T = results['HMP3_T']

ibdmdb = results['ibdmdb']
ibdmdb_T = results['ibdmdb_T']

ibdmdb2 = results['ibdmdb2']
ibdmdb2_T = results['ibdmdb2_T']

Lewis = results['Lewis']
Lewis_T = results['Lewis_T']

MetaHit = results['MetaHit']
MetaHit_T = results['MetaHit_T']

MetaHit2 = results['MetaHit2']
#MetaHit2_T = results['MetaHit2_T']

disease_counts = results['disease_counts']

# 35 Metaphlan control samples are missing
ID_Eran_Elinav_C = Eran_Elinav.loc[Eran_Elinav['Disease'] == 'nonIBD', 'index'].values
Eran_Elinav_C = Eran_Elinav_T.loc[Eran_Elinav_T.index.intersection(ID_Eran_Elinav_C)]

# ID_Gevers_C = Gevers.loc[Gevers["Diagnosis"] == "nonIBD", ...].values - is missing

ID_HMP1_C = HMP1.loc[HMP1["study design"] == 'nonIBD', 'index'].values
HMP1_C = HMP1_T.loc[HMP1_T.index.intersection(ID_HMP1_C)]

ID_HMP2_C = HMP2.loc[HMP1["study design"] == 'nonIBD', 'index'].values
HMP2_C = HMP2_T.loc[HMP2_T.index.intersection(ID_HMP2_C)]


ID_HMP2_pilot_C = HMP2_pilot.loc[HMP2_pilot["diagnosis"] == 'nonIBD', 'index'].values
HMP2_pilot_C = HMP2_pilot_T.loc[HMP2_pilot_T.index.intersection(ID_HMP2_pilot_C)]

ID_HMP3_C = HMP3.loc[HMP3["isolation_source"] == 'GUT', 'index'].values
HMP3_C = HMP3_T.loc[HMP3_T.index.intersection(ID_HMP3_C)]

# 2 Metaphlan control samples are missing
ID_ibdmdb_C = ibdmdb.loc[ibdmdb["diagnosis"] == 'nonIBD', 'index'].values
ibdmdb_C = ibdmdb_T.loc[ibdmdb_T.index.intersection(ID_ibdmdb_C)]

# 1 Metaphlan control samples are missing
ID_ibdmdb2_C = ibdmdb2.loc[ibdmdb2["diagnosis"] == 'nonIBD', 'index'].values
ibdmdb2_C = ibdmdb2_T.loc[ibdmdb2_T.index.intersection(ID_ibdmdb2_C)]

ID_Lewis_C = Lewis.loc[Lewis["diagnosis"] == 'nonIBD', 'index'].values
Lewis_C = Lewis_T.loc[Lewis_T.index.intersection(ID_Lewis_C)]

ID_MetaHit_C = MetaHit.loc[MetaHit["diagnosis"] == 'nonIBD', 'index'].values
MetaHit_C = MetaHit_T.loc[MetaHit_T.index.intersection(ID_MetaHit_C)]

# ID_MetaHit2_C = MetaHit2.loc[MetaHit2["diagnosis"] == "nonIBD", ...].values - is missing

# List of DataFrame names
df_list = [("Eran Elinav" ,Eran_Elinav_C), ("ibdmdb", ibdmdb_C), ("ibdmdb2", ibdmdb2_C), ("HMP1", HMP1_C), ("HMP2",HMP2_C),
          ("HMP2 pilot", HMP2_pilot_C),("HMP3" , HMP3_C), ("Lewis", Lewis_C),("MetaHit", MetaHit_C) ]


def combine_dataframes(df1, df_name1, df2, df_name2):
    df1 = df1.copy()
    df2 = df2.copy()

    df1['SOURCE'] = df_name1
    df2['SOURCE'] = df_name2

    combined_df = pd.concat([df1, df2], axis=0)
    return combined_df


def preprocess_dataframe(df, TAXA_level):
    SOURCE = df['SOURCE'].copy()
    df = df.drop(columns=['SOURCE']).apply(pd.to_numeric, errors='coerce')
    df = df.where(df >= 0.1)
    df = df.fillna(0)
    df = df.loc[:, (df != 0).any(axis=0)]
    df = df[df.sum(axis=1) != 0]

    if TAXA_level == "g":
        df = df.loc[:, df.columns.str.startswith('g__')]
    elif TAXA_level == "s":
        df = df.loc[:, df.columns.str.startswith('s__')]
    df['SOURCE'] = SOURCE
    return df


def run_ANCOM(df, SOURCE, alpha=0.05, tau=0.02, theta=0.1, multiple_comparisons_correction='holm-bonferroni'):
    SOURCE = pd.Series(SOURCE, index=df.index)  # Convert SOURCE to a pandas Series
    replaced_zero = multiplicative_replacement(df.values)
    count_data_replaced_df = pd.DataFrame(replaced_zero, index=df.index, columns=df.columns)
    result = ancom(
        count_data_replaced_df,
        SOURCE,
        alpha=alpha,
        tau=tau,
        theta=theta,
        multiple_comparisons_correction=multiple_comparisons_correction)
    return result[0][result[0]['Reject null hypothesis'] == True][['W']]


def shuffle_labels(labels):
    unique_labels = np.unique(labels)
    shuffled_labels = []
    for label in unique_labels:
        count = np.sum(labels == label)
        shuffled_labels.extend([label] * count)
    np.random.shuffle(shuffled_labels)
    return np.array(shuffled_labels)


def analyze_dataframe(df, cutoff=0):
    SOURCE = df['SOURCE']
    df = df.drop(columns=['SOURCE'])

    # Compute the cutoff using shuffled labels
    random_max_values = []
    for i in range(2):  # Number of iterations to create cutoff
        random_result = run_ANCOM(df, shuffle_labels(SOURCE.values))
        random_max_values.append(random_result['W'].max())

    cutoff = np.nanmax(random_max_values)  # Use nanmax to ignore NaN values
    if np.isnan(cutoff):  # If the cutoff is NaN, replace with default value
        cutoff = 0

    print(f"Cutoff value: {cutoff}")  # Print the computed cutoff

    # Compute the original results using actual labels
    original_results = run_ANCOM(df, SOURCE.values)
    display(original_results)

    # Filter the taxa based on the cutoff
    significant_features = [taxa for taxa, W_value in original_results['W'].items() if W_value > cutoff]
    return significant_features


def analyze_combined_df_list(df_list, TAXA_level):
    all_significant_features = {}
    combinations = itertools.combinations(df_list, 2)

    for (df_name1, df1), (df_name2, df2) in combinations:
        combined_df = combine_dataframes(df1, df_name1, df2, df_name2)
        df_clean = preprocess_dataframe(combined_df, TAXA_level)

        significant_features = analyze_dataframe(df_clean)
        combined_name = f"{df_name1} vs {df_name2}"
        all_significant_features[combined_name] = significant_features
    return all_significant_features


def main():

    result = analyze_combined_df_list(df_list, "g")
    print(result)

if __name__ == '__main__':
    main()