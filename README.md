# neoepitope-nullomers
Workflow to produce nullomers corresponding to neoepitopes from various sources. These scripts also create the alignment file for gene-fusion reads and the nullomer to neoepitope map files

## Prerequisites
- Python (Biopython, numpy)
- R (dplyr, data.table, stringr)
- bedtools
- samtools
- genome fasta and cds fasta, and protein fasta (https://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/)
- TSNAdb_frequent_neoantigen_TCGA_4.0_adj.txt , TSNAdb_frequent_neoantigen_ICGC_4.0_adj.txt
- IEDB_epitope_full_v3.csv, IEDB_mhc_ligand_full.csv, IEDB_tcell_full_v3.csv, IEDB_bcell_full_v3.csv from IEDB database Downloads page
- ChimerKB4.xlsx from ChimerDB. Keep rows with Pub==TRUE and export as a tsv file to ChimerKB_n_Pub.tab
- TE_neoantigens.tsv and TE_chimeric_2297.tsv from Shah et al. 2023 repository

## Create the files required to run cfRNA-neoepitopes pipeline
 
```
createEpitopeDB-nullomers.sh
createFusion-nullomers.sh
createTE-nullomers.sh
```
