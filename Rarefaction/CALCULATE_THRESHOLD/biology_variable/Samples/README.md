# Rarefaction experiments - for biological variation

Process Overview:
   - The first DataFrame is a summary of samples obtained from the European Nucleotide Archive (ENA) at https://www.ebi.ac.uk/ena/browser/view/PRJEB14847. This DataFrame contains detailed information about various samples.
   - The second DataFrame is sourced from the supplementary materials of the mentioned article. This DataFrame provides additional context and description for each sample.
Purpose:
- The script aims to filter and download 41 unique samples based on maximum base count from a total of 21 protocols, each having two stool samples (A & B). It addresses the scenario where the total count is 41 instead of 42 due to a mismatch in one of the samples.

This directory contains the code+results for the rarefaction experiments:
1. Determine what is the correct abunadance threshold for each sequencing depth
2. Determine the change in beta diversity measures based on sequencing depth
3. Determine whether the variation is due to the production processes

**Main** 
The goal of this analysis was to make a clear distinction between the biological variation seen in these samples and the variation noted in our previous simulations that focused on sequence depth. This comparative approach was essential to detect the specific effects of the sample processing, DNA extraction and sequencing procedures apart from the biological variability inherent in the samples. By carefully selecting demographically similar samples, we aimed to reduce the influence of external variables, thus highlighting the true essence of biodiversity.

## Determining abundance thresholds for sequencing depths
### Datasets
**Note** The location Path in severe for all sample in </data1/Oded.Sabah/metaanalysis/biology_variable/Sample_output_bio/PRJEB14847/raw.d>

| Protocol | Sample Accession A | Read Count A | Base Count A | Sample Accession B | Read Count B | Base Count B |
|----------|----------------|-----------|------------|--------------------|--------------|--------------|
| 1        | SAMEA4347635   | 39309536  | 7614194289 | SAMEA4347685       | 43212228     | 8466741542   |
| 2        | SAMEA4347597   | 27318952  | 5354834157 | SAMEA4347592       | 26089883     | 5112772919   |
| 3        | SAMEA4347684   | 31324956  | 6182125911 | SAMEA4347460       | 30403658     | 5921472407   |
| 4        | SAMEA4347581   | 28132202  | 5538422196 | SAMEA4347695       | 33397174     | 6559217672   |
| 5        | SAMEA4347714   | 25948395  | 5022195922 | SAMEA4347692       | 34385115     | 6674743642   |
| 6        | SAMEA4347731   | 28684940  | 5583949512 | SAMEA4347728       | 24901350     | 4939716105   |
| 7        | SAMEA4347601   | 25951297  | 5050425654 | SAMEA4347702       | 28086705     | 5519628633   |
| 8        | SAMEA4347642   | 30011046  | 5779560925 | SAMEA4347654       | 28421572     | 5487704180   |
| 9        | SAMEA4347531   | 29305807  | 5720660469 | SAMEA4347720       | 26291385     | 5153239908   |
| 10       | SAMEA4347723   | 26545834  | 5164753251 | SAMEA4347588       | 24809708     | 4862658854   |
| 11       | SAMEA4347638   | 25047905  | 4845187711 | SAMEA4347678       | 28469684     | 5530213840   |
| 12       | SAMEA4347667   | 27491706  | 5405034177 | SAMEA4347456       | 28782255     | 5659117369   |
| 13       | SAMEA4347533   | 27557387  | 5471358247 | SAMEA4347516       | 28446772     | 5627044342   |
| 14       | SAMEA4347487   | 27632767  | 5351278803 | SAMEA4347480       | 28909611     | 5730326757   |
| 15       | SAMEA4347691   | 26103291  | 5071632686 | SAMEA4347647       | 27746630     | 5380775820   |
| 16       | SAMEA4347530   | 30914915  | 6128535817 | SAMEA4347495       | 28472646     | 5515611139   |
| 17       | SAMEA4347505   | 26912915  | 5309729928 | SAMEA4347721       | 29908015     | 5929156513   |
| 18       | SAMEA4347677   | 28606190  | 5599201768 | SAMEA4347734       | 26028473     | 5152487870   |
| 19       | SAMEA4347625   | 30105212  | 5897286823 | SAMEA4347631       | 26742964     | 5253870213   |
| 20       | SAMEA4347545   | 28729173  | 5624021017 | SAMEA4347523       | 50624665     | 9961003517   |
| 21       |                |           |            | SAMEA4347546       | 42449701     | 8295658430   |


