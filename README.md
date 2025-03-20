
# Rarefaction Analysis for Microbiome Data

## Overview
This project is part of a master's thesis investigating taxa associated with Crohn's disease/Ulcerative colitis in human microbiome samples. The code provides comprehensive tools for rarefaction analysis of metagenomic datasets to determine optimal sequencing depth, taxonomic stability, and abundance thresholds.

## Research Question
The research aims to identify a consistent set of taxa that are associated with Crohn's disease/Ulcerative colitis compared to healthy humans, while determining best practices for metagenomic analysis and taxonomic profiling.

## Components
The pipeline consists of several Python modules working together:

| Module | Function |
|--------|----------|
| `main_cal_thresh_taxa.py` | The main script that orchestrates execution of all components |
| `fastq_selector.py` | Subsamples FASTQ files at different depths (e.g., 1M, 2M, 5M reads) |
| `fastq_splitter.py` | Splits FASTQ files into multiple groups for parallel analysis |
| `Stat_from_select.py` | Analyzes taxonomic profiles from depth-based subsampling |
| `Stat_from_split.py` | Analyzes taxonomic profiles from group-based splitting |

## Data Flow

```
FASTQ Input File
      ↓
 ┌────────────┐
 │ Read Stats │
 └────────────┘
      ↓
┌─────────────────────┐   ┌─────────────────────┐
│  fastq_selector.py  │   │  fastq_splitter.py  │
└─────────────────────┘   └─────────────────────┘
      ↓                          ↓
┌─────────────────────┐   ┌─────────────────────┐
│   MetaPhlAn runs    │   │   MetaPhlAn runs    │
└─────────────────────┘   └─────────────────────┘
      ↓                          ↓
┌─────────────────────┐   ┌─────────────────────┐
│ Stat_from_select.py │   │ Stat_from_split.py  │
└─────────────────────┘   └─────────────────────┘
      ↓                          ↓
    Results                    Results
```

## Usage

```bash
python3 main_cal_thresh_taxa.py <Path to FASTQ file> <Path to Read Stats file> <Output directory> [--thresholds list_of_thresholds]
```

### Required Arguments:
- `<Path to FASTQ file>`: Path to the input FASTQ file (gzipped)
- `<Path to Read Stats file>`: Path to the file containing read statistics
- `<Output directory>`: Directory where all results will be saved

### Optional Arguments:
- `--thresholds`: List of abundance thresholds for filtering taxa (default: [0, 0.001, 0.01, 0.1, 1])

### Example:
```bash
python3 main_cal_thresh_taxa.py /data2/Human/HMP.phase1/Stool/raw.d/SAMN00040286/SRS019068.denovo_duplicates_marked.trimmed.1.fastq.gz /data2/Human/HMP.phase1/Stool/raw.d/Stool.read-stats.txt /data1/Oded.Sabah/metaanalysis/SAMN00040286
```

## Pipeline Details

### 1. Input Processing
The pipeline starts by calculating the total number of reads in the sample from the read stats file. This is used to determine the subsampling strategy.

### 2. Depth Analysis (`fastq_selector.py`)
This module selects different numbers of reads from the input FASTQ at specific increments:

| Read Depth | Description |
|------------|-------------|
| 1M - 9M    | Every 1M reads |
| 10M - 18M  | Every 2M reads |
| 20M - 45M  | Every 5M reads |
| 50M - 90M  | Every 10M reads |
| 100M - 1000M | Every 100M reads |

For each subsample:
- Reads are selected systematically by using a "jump" value
- A new FASTQ file is created
- MetaPhlAn is run on the subsample
- Output taxonomic profiles are saved

### 3. Split Analysis (`fastq_splitter.py`)
This module takes a different approach by:
- Incrementing the read count by 1M reads in each iteration
- For each count, splitting the data into 10 different groups
- Running MetaPhlAn on each group separately
- This allows assessment of variability within the same sample at different depths

### 4. Statistical Analysis

#### Depth-based Analysis (`Stat_from_select.py`)
For each abundance threshold (0%, 0.001%, 0.01%, 0.1%, 1%):
- Calculates alpha diversity metrics:
  - Richness (number of observed species)
  - Shannon diversity index
  - Simpson's index and evenness
- Computes Jaccard similarity between each subsample and the full sample (ground truth)
- Calculates Bray-Curtis dissimilarity, L1, and L2 distances

#### Split-based Analysis (`Stat_from_split.py`)
For each abundance threshold:
- Compares taxonomic profiles between different split groups
- Calculates Jaccard similarity between group pairs
- Computes Bray-Curtis dissimilarity, L1, and L2 distances
- Evaluates consistency across splits at the same read depth

## Output Files

### Directory Structure
```
<Output Directory>/
├── mp4.depth.<sample_id>/
│   ├── <sample_id>.<depth>M.fastq.gz          # Subsampled FASTQ files
│   ├── <sample_id>.<depth>M.txt               # MetaPhlAn output files
│   └── <sample_id>.<depth>M.txt.stdout        # MetaPhlAn log files
├── mp4.split.<sample_id>/
│   ├── <sample_id>_group_<N>-<depth>M.fastq.gz  # Split FASTQ files
│   ├── <sample_id>_group_<N>-<depth>M.txt       # MetaPhlAn output for splits
│   └── <sample_id>_group_<N>-<depth>M.txt.stdout # MetaPhlAn logs for splits
├── <sample_id>_thresh_<threshold>.csv         # Results from depth analysis
└── <sample_id>_split_<threshold>.csv          # Results from split analysis
```

### Result Files
For each threshold value, CSV files are generated containing:

#### Depth Analysis Results:
- Column name (read depth)
- Number of species detected
- Jaccard similarity to ground truth
- Richness, Diversity, and Evenness metrics
- Bray-Curtis dissimilarity
- L1 and L2 distance scores

#### Split Analysis Results:
- Pair of columns compared
- Mean number of species
- Jaccard similarity between pairs
- Bray-Curtis dissimilarity
- L1 and L2 distance scores

## Technical Details

### Depth Selection Strategy
The `selected_reads` method in `FastqProcessor` determines which sequencing depths to analyze:
```python
n_serial = (list(range(1, 10)) + list(range(10, 19, 2)) +
            list(range(20, 49, 5)) + list(range(50, 99, 10)) +
            list(range(100, 1001, 100)))
n_serial = [i * 1e6 for i in n_serial]
```

### Subsampling Mechanism
Reads are selected using a "jump" value. For example, to select 1M reads from a 50M sample:
- Jump = 50 (select every 50th read)
- This ensures even distribution across the file

### Taxonomy Filtering
The pipeline processes taxonomic data at multiple threshold levels (default: 0%, 0.001%, 0.01%, 0.1%, 1%), allowing examination of how abundance thresholds affect taxonomic profiles.

## Requirements

### Software Requirements
- Python 3.6+
- pandas
- scikit-bio
- gzip
- MetaPhlAn 4

### Input Requirements
- Gzipped FASTQ file
- Read statistics file with tab-delimited format containing sample IDs and read counts

## Research Applications
This tool addresses several key questions in microbiome research:
1. What is the minimum sequencing depth required for reliable taxonomic profiling?
2. How do abundance thresholds affect the stability of taxonomic profiles?
3. How confident can we be in the taxa identified at different sequencing depths?
4. What is the variation between technical replicates at the same sequencing depth?

The insights gained from this analysis help establish benchmarks for future metagenomic studies, particularly those focused on identifying disease-associated taxa.
