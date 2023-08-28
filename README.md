# neoepitope-nullomers
Workflow to produce nullomers corresponding to neoepitopes from various sources. These scripts also create the alignment file for gene-fusion reads and the nullomer to neoepitope map files

## Prerequisites
Python (Biopython, numpy)
R (dplyr, data.table, stringr)
bedtools
samtools
genome fasta and cds fasta, and protein fasta (https://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/)
TSNAdb_frequent_neoantigen_TCGA_4.0_adj.txt , TSNAdb_frequent_neoantigen_ICGC_4.0_adj.txt

## Create all these files
 
```
createEpitopeDB-nullomers.sh
createFusion-nullomers.sh
createTE-nullomers.sh
```
