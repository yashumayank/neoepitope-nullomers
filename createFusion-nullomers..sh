#!/bin/bash
#BSUB -J epiDB-Null
#BSUB -o test-%J.out
#BSUB -e test-%J.err
#BSUB -q bigmem
#BSUB -n 1

#eval "$(/PHShome/mm1169/miniconda3/bin/conda shell.bash hook)"
#conda activate scanpy
cd /data/hemberg/nullomers/IEDB/epitope_DBs/fusion
#--extract Human cds 16mers
#awk '{if($0~/^>.*/){ln=length(x);for(i=1;i+15<=ln;i++){a[substr(x,i,16)]=1};x=""}else{x=x $1}}END{for(j in a){print j}}' ../Homo_sapiens.GRCh38.cds.all.fa > Homo_sapiens_cds_16mers.tab

#1342 have valid junctions out of 2243 fusions that are part of both ChimerKB & Pub
#9(17)-15(577) cds nullomers per junction
#0(16)-15(157) genome nullomers per junction 2(28), 3(35)

#create protein kmer to get the null peptides
python3 ../createProteinKmers.py  ../../Homo_sapiens.GRCh38.pep.all.fa.gz > ../Homo_sapiens_protein_9mers.tab

zcat /data/hemberg/shared_resources/genomes/human/hg19.fa.gz|awk '{if($0~/^>.*/){ln=length(x);for(i=1;i+15<=ln;i++){a[substr(x,i,16)]=1};x=""}else{x=x $1}}END{for(j in a){print j}}' - > ../Homo_sapiens_genome_16mers.tab

python3 extractChimerDbNullomers.py ChimerKB_n_Pub.tab /data/hemberg/shared_resources/genomes/human/hg19.fa.gz /data/hemberg/shared_resources/genomes/human/hg19.ensGene.gtf ChimerKB_n_Pub
#Fusion junctions ouside CDS: 680
python3 ChimerKB_5_3_seqs2nullomers.py /data/hemberg/nullomers/IEDB/epitope_DBs/Homo_sapiens_cds_16mers.tab ../../Homo_sapiens.GRCh38.pep.all.fa.gz ChimerKB_n_Pub

awk '{split($1,u,";");for(i=1;i<length(u);i++)print toupper(u[i])}' ChimerKB_n_Pub_junction_nullomers.fasta > ChimerKB_nullomers.tsv
#revComp transcriptome_nullomers to search on the read 2
awk 'BEGIN{c["A"]="T";c["T"]="A";c["G"]="C";c["C"]="G";}{y="";for(j=16;j>=1;j--){y=y c[substr($1,j,1)]};print y;}' ChimerKB_nullomers.tsv > ChimerKB_nullomers_revComp.tsv
