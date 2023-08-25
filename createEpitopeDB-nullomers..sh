#!/bin/bash
#BSUB -J epiDB-Null
#BSUB -o test-%J.out
#BSUB -e test-%J.err
#BSUB -q bigmem
#BSUB -n 1

eval "$(/PHShome/mm1169/miniconda3/bin/conda shell.bash hook)"
conda activate scanpy
cd /data/hemberg/nullomers/IEDB/epitope_DBs/
#--extract Human cds 16mers
#awk '{if($0~/^>.*/){ln=length(x);for(i=1;i+15<=ln;i++){a[substr(x,i,16)]=1};x=""}else{x=x $1}}END{for(j in a){print j}}' ../Homo_sapiens.GRCh38.cds.all.fa > Homo_sapiens_cds_16mers.tab

#--extract extract nullomers associated to neoepitopes in IEDB and TSNAdb
python IEDB_TSNAdb2nullomer.py TSNAdb_frequent_ICGC_per_ENST3.tab Homo_sapiens_cds_16mers.tab ../Homo_sapiens.GRCh38.cds.all.fa TSNAdb_ICGC
python IEDB_TSNAdb2nullomer.py TSNAdb_frequent_TCGA_per_ENST3.tab Homo_sapiens_cds_16mers.tab ../Homo_sapiens.GRCh38.cds.all.fa TSNAdb_TCGA
python IEDB_TSNAdb2nullomer.py IEDB_neoepitopes_per_ENST2.tab Homo_sapiens_cds_16mers.tab ../Homo_sapiens.GRCh38.cds.all.fa IEDB

#---Count the WT-neoepitope pairs extracted from each database
echo "TSNAdb_ICGC"
awk 'NF==4' TSNAdb_ICGC_neoepitopes-nullomers.tsv |cut -f3 |sort|uniq -c|sort|wc -l
echo "TSNAdb_TCGA"
awk 'NF==4' TSNAdb_TCGA_neoepitopes-nullomers.tsv |cut -f3 |sort|uniq -c|sort|wc -l
echo "IEDB"
awk 'NF==4' IEDB_neoepitopes-nullomers.tsv |cut -f3 |sort|uniq -c|sort|wc -l

#--merge the nullomers associated to all the above neoepitopes
awk -F "\t" '{split($1,u,";");for(i in u){a[u[i]]=1}}END{for(i in a)print i}' IEDB_neoepitopes-nullomers.tsv TSNAdb_TCGA_neoepitopes-nullomers.tsv TSNAdb_ICGC_neoepitopes-nullomers.tsv > epitopeDB_nullomers.all.tsv

#--Calculate genome frequency of each nullomer
awk -F "\t" '(NR==FNR){a[$1]=1;next}{if($1~/^>.*/){if(seqid!=""){for(i=1;i<=(length(seq)-15);i++){xtr=substr(seq,i,16);if(xtr in a){a[xtr]++}}};seqid=$1;seq=""}else{seq=seq""$1}}END{for(j=1;j<=(length(seq)-15);j++){xtr=substr(seq,j,16);if(xtr in a){a[xtr]++}};for(k in a){print k"\t"a[k]-1}}' epitopeDB_nullomers.all.tsv /data/hemberg/shared_resources/genomes/human/GRCh38.p13.107.genome.fa > epitopeDB_nullomers_freq.tsv

#---take 20 least frequent per neoepitope with genome-frequency < 100
awk -F "\t" 'BEGIN{PROCINFO["sorted_in"] = "@ind_num_asc"}{if(NR==FNR){f[$1]=$2}else{split($1,u,";");delete(v);delete(vX);for(i=1;i<=length(u);i++){v[i] = f[u[i]]};n = asorti(v,vX,"@val_num_asc");c=1;nlist="";for(x=1;x<=100;x++){while(v[vX[c]]<x && c<=n){nlist=nlist";"u[vX[c]];c++};if(c>20 || c>n){break}};print substr(nlist,2)"\t"$2"\t"$3"\t"$4"\t"$5}}' epitopeDB_nullomers_freq.tsv TSNAdb_ICGC_neoepitopes-nullomers.tsv > TSNAdb_ICGC_neoepitopes-nullomersTop20.tsv
awk -F "\t" 'BEGIN{PROCINFO["sorted_in"] = "@ind_num_asc"}{if(NR==FNR){f[$1]=$2}else{split($1,u,";");delete(v);delete(vX);for(i=1;i<=length(u);i++){v[i] = f[u[i]]};n = asorti(v,vX,"@val_num_asc");c=1;nlist="";for(x=1;x<=100;x++){while(v[vX[c]]<x && c<=n){nlist=nlist";"u[vX[c]];c++};if(c>20 || c>n){break}};print substr(nlist,2)"\t"$2"\t"$3"\t"$4"\t"$5}}' epitopeDB_nullomers_freq.tsv TSNAdb_TCGA_neoepitopes-nullomers.tsv > TSNAdb_TCGA_neoepitopes-nullomersTop20.tsv
awk -F "\t" 'BEGIN{PROCINFO["sorted_in"] = "@ind_num_asc"}{if(NR==FNR){f[$1]=$2}else{split($1,u,";");delete(v);delete(vX);for(i=1;i<=length(u);i++){v[i] = f[u[i]]};n = asorti(v,vX,"@val_num_asc");c=1;nlist="";for(x=1;x<=100;x++){while(v[vX[c]]<x && c<=n){nlist=nlist";"u[vX[c]];c++};if(c>20 || c>n){break}};print substr(nlist,2)"\t"$2"\t"$3"\t"$4"\t"$5}}' epitopeDB_nullomers_freq.tsv IEDB_neoepitopes-nullomers.tsv > IEDB_neoepitopes-nullomersTop20.tsv

#--merge the filtered nullomers associated to IEDB and TSNAdb neoepitopes
awk -F "\t" '($1!="" || NF==4){split($1,u,";");for(i in u){a[u[i]]=1}}END{for(i in a)print i}' IEDB_neoepitopes-nullomersTop20.tsv TSNAdb_TCGA_neoepitopes-nullomersTop20.tsv TSNAdb_ICGC_neoepitopes-nullomersTop20.tsv > epitopeDB_nullomers.tsv
#--revComp cds_nullomers to search on the read 2 of the RNAseq data
awk 'BEGIN{c["A"]="T";c["T"]="A";c["G"]="C";c["C"]="G";}{y="";for(j=16;j>=1;j--){y=y c[substr($1,j,1)]};print y;}' epitopeDB_nullomers.tsv > epitopeDB_nullomers_revComp.tsv

sort TSNAdb_ICGC_neoepitopes-nullomersTop20.tsv|uniq >TSNAdb_ICGC_neoepitopes-nullomers.filtered.tsv
sort TSNAdb_TCGA_neoepitopes-nullomersTop20.tsv|uniq >TSNAdb_TCGA_neoepitopes-nullomers.filtered.tsv
sort IEDB_neoepitopes-nullomersTop20.tsv|uniq > IEDB_neoepitopes-nullomers.filtered.tsv

