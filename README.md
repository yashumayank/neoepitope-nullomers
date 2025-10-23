# neoepitope-nullomers
Workflow to extract known human neoepitopes from various sources and produce nullomers corresponding to those neoepitopes. For each known neoepitope, their affinity of common HLA alleles in estimated using netMHCpan and there prevelance in the healthy human germline is estimated using GnomAD database. These scripts also create the alignment files and the nullomer-neoepitope map files that are used by the FastNeo (https://academic.oup.com/bioinformatics/article/41/5/btaf138/8124074)

## Prerequisites
- Python (Biopython, numpy)
- R (dplyr, data.table, stringr)
- bedtools
- samtools
- genome fasta and cds fasta, and protein fasta (https://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/)
- TSNAdb_frequent_neoantigen_TCGA_4.0_adj.txt, TSNAdb_frequent_neoantigen_ICGC_4.0_adj.txt from TSNAdb (http://biopharm.zju.edu.cn/tsnadb/download/)
- epitope_full_v3.csv, mhc_ligand_full.csv, tcell_full_v3.csv, bcell_full_v3.csv from IEDB database (downloaded on Sept-19-2022)
- ChimerKB4.xlsx from ChimerDB. Keep rows common with 'Pub' and export as a tsv file to ChimerKB_n_Pub.tab  (downloaded on May-19-2023)
- Genome fasta and gtf file from hg19 (https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/)
- TE_neoantigens.tsv and TE_chimeric_2297.tsv from Shah et al. 2023 repository
- protein ID to transcript id mapiing files from uniprot and gencode HUMAN_9606_idmapping_selected.tab, gencode.v40.metadata.TrEMBL, gencode.v40.metadata.SwissProt  (downloaded on June-27-2022)

## Pipelines to create neoepitope associated nullomers 
The following pipelines create nullomers associated to the known neoepitopes or gene fusions. A nullomer mapping file is created by each of the scripts. The neoepitope-nullomer mapping file also includes the differential HLA affinity of the neoepitopes and their presence among the germline variants. Recently created mapping files can be found in the FastNeo repository.
 
```
createEpitopeDB-nullomers.sh
createFusion-nullomers.sh
```
