# Rarefaction Analysis Software

## Overview

This software package provides tools for conducting rarefaction experiments on metagenomic sequencing data. The primary goals are:

1. Determine the optimal abundance threshold for taxa at different sequencing depths
2. Measure changes in beta diversity metrics based on sequencing depth

## Architecture

The software consists of several Python modules that work together to process FASTQ files, run Metaphlan analysis at different read depths, and calculate diversity statistics.

```
├── CALCULATE_THRESHOLD/
│   ├── fastq_selector.py     # Selects subsets of reads at different depths
│   ├── fastq_splitter.py     # Splits FASTQ files into multiple groups
│   ├── Stat_from_select.py   # Analyzes data from selected reads
│   ├── Stat_from_split.py    # Analyzes data from split reads
│   └── main_cal_thresh_taxa.py  # Main execution script
```

## Workflow

The analysis pipeline follows these steps:

1. **Read Selection/Splitting**: The software either selects subsets of reads from a FASTQ file at different depths (`fastq_selector.py`) or splits the file into multiple groups (`fastq_splitter.py`).

2. **Metaphlan Analysis**: For each subset or group, Metaphlan is run to identify and quantify microbial taxa.

3. **Statistical Analysis**: The resulting data is analyzed to calculate:
   - Alpha diversity metrics (richness, Shannon diversity, Simpson's evenness)
   - Beta diversity metrics (Jaccard similarity, Bray-Curtis dissimilarity)
   - L1 and L2 distance scores

4. **Results**: The analysis generates CSV files with the calculated statistics for different abundance thresholds.

## Components

### FastqProcessor

The `FastqProcessor` class in `fastq_selector.py` selects subsets of reads from a FASTQ file at different depths.

Key functionalities:
- Calculates appropriate sample sizes based on the input FASTQ file
- Creates subsets at various read depths (1M, 2M, 3M, ..., 1000M reads)
- Runs Metaphlan on each subset to identify microbial taxa

### FastqSplitter

The `FastqSplitter` class in `fastq_splitter.py` splits a FASTQ file into multiple groups for comparative analysis.

Key functionalities:
- Divides reads into 10 equal groups at different depths
- Processes each group separately through Metaphlan
- Enables comparison between groups at the same depth

### MetaphlanAnalysis

The `MetaphlanAnalysis` class in `Stat_from_select.py` analyzes the output from `FastqProcessor`.

Key functionalities:
- Computes alpha diversity scores (richness, Shannon diversity, Simpson's evenness)
- Calculates Jaccard similarity scores against the ground truth (sample with maximum number of readings)
- Measures Bray-Curtis dissimilarity and L1/L2 distances

### MetaphlanAnalysis_from_split

The `MetaphlanAnalysis_from_split` class in `Stat_from_split.py` analyzes the output from `FastqSplitter`.

Key functionalities:
- Computes Jaccard similarity scores between paired samples
- Calculates Bray-Curtis dissimilarity and L1/L2 distances between sample pairs
- Sorts and formats results for better interpretation

## Usage

The main script `main_cal_thresh_taxa.py` coordinates the entire workflow and can be executed as follows:

```bash
python3 main_cal_thresh_taxa.py <path_to_fastq_file> <path_to_read_stats> <output_directory> [--thresholds <threshold_list>]
```

Example:
```bash
python3 main_cal_thresh_taxa.py /data2/Human/HMP.phase1/Stool/raw.d/SAMN00040286/SRS019068.denovo_duplicates_marked.trimmed.1.fastq.gz /data2/Human/HMP.phase1/Stool/raw.d/Stool.read-stats.txt /data1/Oded.Sabah/metaanalysis/SAMN00040286
```

### Parameters:

- `fastq_file`: Path to the input FASTQ file
- `read_stats`: Path to the file containing read statistics
- `output_dir`: Directory where results will be stored
- `--thresholds`: Optional list of abundance thresholds (default: [0, 0.001, 0.01, 0.1, 1])

## Input Data

The software has been tested with high-coverage metagenomic samples, mainly from human microbiome projects:

| Dataset | Samples | Sample Size (reads) |
|---------|---------|---------------------|
| HMP Phase 1 Stool | 3 | ~146-149M |
| HMP Phase 2 Stool | 3 | ~112-131M |
| PRJEB50555 | 2 | ~189-258M |

## Output Files and Results

The software generates several output files in the specified output directory:

### 1. Depth Analysis Results

File: `<id_sample>_thresh_<threshold>.csv`

This file contains the analysis results from the `FastqProcessor` and `MetaphlanAnalysis` components, with metrics for each sequencing depth against the ground truth (maximum depth sample).

Example content:

| Col_name | Num_of_species | Ground_truth | Richness | Diversity | Evenness | Dissimilarity | l1_score | l2_score |
|----------|----------------|--------------|----------|-----------|----------|---------------|----------|----------|
| sample.1M | 42 | 0.65 | 42.0 | 3.21 | 0.78 | 0.24 | 0.35 | 0.12 |
| sample.2M | 58 | 0.72 | 58.0 | 3.45 | 0.82 | 0.18 | 0.28 | 0.09 |
| sample.5M | 81 | 0.85 | 81.0 | 3.68 | 0.85 | 0.12 | 0.19 | 0.06 |
| sample.10M | 93 | 0.91 | 93.0 | 3.74 | 0.88 | 0.07 | 0.12 | 0.04 |
| sample.20M | 98 | 0.95 | 98.0 | 3.78 | 0.89 | 0.04 | 0.08 | 0.03 |
| sample.50M | 102 | 0.97 | 102.0 | 3.80 | 0.90 | 0.02 | 0.05 | 0.02 |
| sample.100M | 105 | 1.00 | 105.0 | 3.81 | 0.91 | 0.00 | 0.00 | 0.00 |

Where:
- `Col_name`: Sample name with read depth (e.g., 1M = 1 million reads)
- `Num_of_species`: Number of detected species at the given depth
- `Ground_truth`: Jaccard similarity score compared to the ground truth sample
- `Richness`, `Diversity`, `Evenness`: Alpha diversity metrics
- `Dissimilarity`: Bray-Curtis dissimilarity from ground truth
- `l1_score`, `l2_score`: L1 and L2 distance metrics

### 2. Split Analysis Results

File: `<id_sample>_split_<threshold>.csv`

This file contains the results from the `FastqSplitter` and `MetaphlanAnalysis_from_split` components, showing comparisons between different groups at the same sequencing depth.

Example content:

| Cols_name | jaccard_scores | m_species | Dissimilarity | l1_score | l2_score |
|-----------|----------------|-----------|---------------|----------|----------|
| group_1-1M VS group_2-1M | 0.625 | 40.5 | 0.285 | 0.412 | 0.152 |
| group_3-1M VS group_4-1M | 0.608 | 39.0 | 0.302 | 0.438 | 0.164 |
| group_1-2M VS group_2-2M | 0.681 | 56.5 | 0.215 | 0.322 | 0.128 |
| group_3-2M VS group_4-2M | 0.695 | 57.0 | 0.198 | 0.305 | 0.118 |
| group_1-5M VS group_2-5M | 0.778 | 80.0 | 0.142 | 0.225 | 0.085 |
| group_1-10M VS group_2-10M | 0.855 | 92.5 | 0.081 | 0.135 | 0.052 |
| group_1-20M VS group_2-20M | 0.904 | 97.0 | 0.052 | 0.089 | 0.032 |

Where:
- `Cols_name`: Comparison between two groups at the same depth
- `jaccard_scores`: Jaccard similarity between the groups
- `m_species`: Average number of species detected in the two groups
- `Dissimilarity`: Bray-Curtis dissimilarity between the groups
- `l1_score`, `l2_score`: L1 and L2 distance metrics

### Visualization and Interpretation

The CSV results can be plotted to visualize important trends:

1. **Rarefaction Curves**: Plotting `Num_of_species` against read depth shows how species richness increases with sequencing depth and eventually plateaus.

2. **Similarity to Ground Truth**: Plotting `Ground_truth` values against read depth illustrates how quickly the taxonomic profile stabilizes.

3. **Dissimilarity Trends**: Plotting `Dissimilarity` against read depth helps determine the minimum sequencing depth required for stable results.

4. **Abundance Threshold Impact**: Comparing results across different threshold values (0, 0.001, 0.01, 0.1, 1) reveals how filtering low-abundance taxa affects stability and diversity metrics.

These analyses help researchers determine:
- Optimal sequencing depth for their metagenomic studies
- Appropriate abundance thresholds for filtering noise
- Minimum sequencing requirements for reliable beta diversity comparisons

## Statistical Measures

### Alpha Diversity Metrics
- **Richness**: Number of observed species (observed OTUs)
- **Shannon Diversity**: Measure of both richness and evenness
- **Simpson's Evenness**: How evenly taxa are distributed in the sample

### Beta Diversity Metrics
- **Jaccard Similarity**: Ratio of the intersection to the union of taxa sets
- **Bray-Curtis Dissimilarity**: Compositional dissimilarity based on abundance
- **L1 Score**: Sum of absolute differences between abundance profiles
- **L2 Score**: Euclidean distance between abundance profiles

## Advanced Features

- **Parallel Processing**: Uses `ProcessPoolExecutor` for concurrent execution of tasks
- **Adaptive Sampling**: Adjusts sampling rates based on input file size
- **Configurable Thresholds**: Allows testing multiple abundance thresholds

## Installation and Setup


### Installation Steps

1. Clone the repository:
   ```bash
   git clone https://github.com/username/rarefaction-analysis.git
   cd rarefaction-analysis
   ```

2. Install required Python packages:
   ```bash
   pip install -r requirements.txt
   ```

   The requirements.txt file contains:
   ```
   pandas>=1.0.0
   scikit-bio>=0.5.6
   numpy>=1.18.0
   ```

   Note: Other dependencies (gzip, concurrent.futures, logging, math, argparse) are part of the Python standard library.

3. Ensure MetaPhlAn is correctly configured:
   - The software expects to find MetaPhlAn at `/data1/software/metaphlan/run-metaphlan.sh`
   - If your installation is different, update the path in `fastq_selector.py` and `fastq_splitter.py`

4. Verify installation:
   ```bash
   python3 CALCULATE_THRESHOLD/main_cal_thresh_taxa.py --help
   ```
   
   You should see the help message with available command-line options.

### Directory Structure

After installation, ensure your directory structure looks like:
```
rarefaction-analysis/
├── CALCULATE_THRESHOLD/
│   ├── fastq_selector.py
│   ├── fastq_splitter.py
│   ├── Stat_from_select.py
│   ├── Stat_from_split.py
│   └── main_cal_thresh_taxa.py
├── requirements.txt
├── README.md
└── REQUESTS.md
```