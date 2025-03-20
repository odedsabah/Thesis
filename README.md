# Comprehensive Analysis of Oded's M.Sc. Thesis

## Thesis Overview

This thesis investigates the relationship between Inflammatory Bowel Disease (IBD) - specifically Crohn's Disease (CD) and Ulcerative Colitis (UC) - and the human gut microbiome. The primary research question is to identify a consistent set of taxa that distinguish between IBD patients and healthy controls across multiple datasets.

## Research Questions and Objectives

The thesis aims to:

1. Determine best practices for comparing phenotypes in microbiome studies
2. Apply these practices to available IBD datasets to identify significant taxa
3. Establish methodological guidelines for future microbiome studies

## Methodological Considerations

### Best Practices Being Explored

| Practice | Key Question | Approach |
|----------|--------------|----------|
| Rarefaction | Should rarefaction be performed for metagenomics (MGX) datasets? | Yes (work already completed) |
| Metadata Effects | How do factors like age, sex, country, and calprotectin influence the microbiome? | Statistical analysis of control groups |
| Dataset Comparability | Are different datasets comparable? | Beta diversity analysis between matched samples |
| Feature Significance | What is a reliable method for determining statistically significant taxa? | Method-agnostic approach with randomization tests |

### Evaluating Metadata Effects on Microbiome

The thesis proposes the following analysis workflow:

1. For each control group in each dataset:
   - Match samples based on features (age, sex, country)
   - Calculate Jaccard and Bray-Curtis distance distributions
   - Compare against random pairs as baseline
2. Visualization approach:
   - Box plots for each feature
   - Histograms for age distribution
   - Bar plots for sex and country distributions
3. Statistical testing:
   - Test significance of differences between feature-matched and random pairs
   - Consider anti-matching to amplify effects
4. Extended analysis:
   - Examine PCoA plots for unexpected clustering
   - Repeat analysis for Crohn's disease group

### Assessing Dataset Comparability

For each pair of datasets (A and B):

1. Match samples based on important features:
   - Cross-dataset pairs (a, b) where a ∈ A and b ∈ B
   - Within-dataset pairs (a₁, a₂) where a₁, a₂ ∈ A
   - Within-dataset pairs (b₁, b₂) where b₁, b₂ ∈ B
2. Draw boxplots of beta diversity distributions
3. Test if between-dataset diversity is significantly different from within-dataset diversity

### Finding Significant Features

The thesis proposes this methodology:

1. Run feature selection method on data with original phenotype labels
2. Repeat 10 times:
   - Randomly split case and control groups in original proportions
   - Run feature selection method
   - Find highest random score and use as cutoff
3. Compare results across multiple feature selection methods

## Data Sources and Datasets

### Potential Datasets for Analysis

![Dataset Relationships](datasets/circle_graph.png)

*Figure 1: Visualization of relationships between the analyzed datasets*

### Metadata to Collect

For each sample, the thesis aims to collect:

| Metadata Category | Description |
|-------------------|-------------|
| Project ID | Dataset identifier |
| Sample ID | Unique sample identifier |
| Subject ID | Subject/patient identifier |
| Time | Sampling timepoint |
| Sequencing Depth | Number of reads and base pairs |
| Health Status | Healthy, CD, UC, or general IBD |
| Demographics | Age, sex, country |
| Clinical Markers | Calprotectin (FCP) |
| Sample Quality | Human DNA percentage |
| Medication | Antibiotics use, steroids use |

## Analysis Pipeline

### High-Level Analysis Plan

1. Literature review of similar meta-analyses
2. Dataset identification and metadata collection
3. Understanding prior approaches to taxa list generation
4. Determining community composition for each sample
5. Assessing dataset comparability
6. Establishing a bioinformatics/statistical procedure
7. Results analysis and conclusions

### Low-Level Analysis Process

#### 1. Community Composition Analysis

**Tool Selection:**
- Using MetaPhlAn4 due to:
  - Popularity and widespread use
  - Most updated reference database
  - Better precision compared to Kraken2/Centrifuge
  - More comprehensive than mOTUs2

**Rarefaction Analysis:**
- Using large samples (e.g., SAMEA110452975 with 299,380,873 reads)
- Creating subsamples at multiple depths (1M to 70M+ reads)
- Evaluating Jaccard index at different abundance thresholds (1%, 0.1%, 0.01%)
- Measuring differences in relative abundance between subsamples and complete samples

**Example Command:**
```bash
conda activate base
/home/itaish/software/processing/run-metaphlan4.sh
```

#### 2. Dataset Comparability

**Naive Approaches:**
1. Finding separating taxa between healthy and IBD groups in each dataset
2. Analyzing differences in community composition (beta-diversity)

**Advanced Approach:**
- Considering confounding parameters (age, sex, health status)

#### 3. Developing Guidelines

The thesis aims to establish guidelines for future studies:
- Minimum number of participants required
- Minimum abundance threshold for meaningful taxa
- Optimal taxonomic levels for observing differences

## Literature Review Focus

The thesis plans to analyze prior meta-analyses with focus on:

| Question | Details to Extract |
|----------|-------------------|
| Data Types | MGX or 16S rRNA sequencing |
| Methodologies | Community profiling approach, taxonomic levels analyzed, statistical/ML methods |
| Comparability | How authors ensured datasets were comparable |
| Analysis Approach | Individual datasets or combined analysis |
| Dataset Weighting | Were datasets weighted by size or equally sampled? |
| Time Series | How were time series datasets handled? |
| Metadata | Were factors like age, sex, nationality, calprotectin considered? |

## Computational Resources

The thesis utilizes server resources with data located in specific directories:
- Raw data: `/data2/Human/[Project]/Stool/raw.d/`
- Sample statistics: `[Project].read-stats.txt`

## Expected Outcomes

1. A consistent set of taxa associated with IBD
2. Methodological guidelines for microbiome comparison studies
3. Insights into required sequencing depth and sample sizes
4. Framework for handling confounding variables in microbiome studies