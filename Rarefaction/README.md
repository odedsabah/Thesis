# Rarefaction experiments

This directory contains the code+results for rhe rarefaction experiments:
1. Determine what is the correct abunadance threshold for each sequencing depth
2. Determine the change in beta diversity measures based on sequencing depth

## Determining abundance thresholds for sequencing depths
### 1. Datasets

We used this sample based on >100M paired-reads and no more than 3 samples from one dataset:

| MetaPhlAn4 Path file                                          | Read Stats                                                | FASTQ File                                                                                                         | n samples | Sample size | ... |
|:--------------------------------------------------------------|:----------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------| :--- | :--- | :--- |
| /data2/Human/HMP.phase1/Stool/metaphlan4.0.6/SAMN00040286.txt | /data2/Human/HMP.phase1/Stool/raw.d/Stool.read-stats.txt  | /data2/Human/HMP.phase1/Stool/raw.d/SAMN00040286/SRS019068.denovo_duplicates_marked.trimmed.1.fastq.gz             | 1 | 149762007 | 14965582301 |
| /data2/Human/HMP.phase1/Stool/metaphlan4.0.6/SAMN00036649.txt | /data2/Human/HMP.phase1/Stool/raw.d/Stool.read-stats.txt  | /data2/Human/HMP.phase1/Stool/raw.d/SAMN00036649/SRS015431.denovo_duplicates_marked.trimmed.1.fastq.gz             | 1 | 149162797 | 14901563545 |
| /data2/Human/HMP.phase1/Stool/metaphlan4.0.6/SAMN00037108.txt | /data2/Human/HMP.phase1/Stool/raw.d/Stool.read-stats.txt  | /data2/Human/HMP.phase1/Stool/raw.d/SAMN00037108/SRS015890.denovo_duplicates_marked.trimmed.1.fastq.gz             | 1 | 146067069 | 14596867636 |
| /data2/Human/HMP.phase2/Stool/metaphlan4.0.6/SAMN00146983.txt | /data2/Human/HMP.phase2/Stool/raw.d/Stool.read-stats.txt  | /data2/Human/HMP.phase2/Stool/raw.d/SAMN00146983/SRS147766.denovo_duplicates_marked.trimmed.1.fastq.gz             | 1 | 131485676 | 13148567600 |
| /data2/Human/HMP.phase2/Stool/metaphlan4.0.6/SAMN00099612.txt | /data2/Human/HMP.phase2/Stool/raw.d/Stool.read-stats.txt  | /data2/Human/HMP.phase2/Stool/raw.d/SAMN00099612/SRS103987.denovo_duplicates_marked.trimmed.1.fastq.gz             | 1 | 116716792 | 11671679200 |
| /data2/Human/HMP.phase2/Stool/metaphlan4.0.6/SAMN00143093.txt | /data2/Human/HMP.phase2/Stool/raw.d/Stool.read-stats.txt  | /data2/Human/HMP.phase2/Stool/raw.d/SAMN00143093/SRS143876.denovo_duplicates_marked.trimmed.1.fastq.gz             | 1 | 112287524 | 11228752400 |
| /data2/Human/PRJEB50555/metaphlan4.0.6/SAMEA110452918.txt     | /data2/Human/PRJEB50555/raw.d/PRJEB50555.read-stats.txt   | /data2/Human/PRJEB50555/raw.d/SAMEA110452918/9835_R1.fastq.gz                                                      | 1 | 258362901 | 38720819401 |
| /data2/Human/PRJEB50555/metaphlan4.0.6/SAMEA110452924.txt     | /data2/Human/PRJEB50555/raw.d/PRJEB50555.read-stats.txt   | /data2/Human/PRJEB50555/raw.d/SAMEA110452924/9883_R1.fastq.gz                                                      | 1 | 189178947 | 28442121817 |

