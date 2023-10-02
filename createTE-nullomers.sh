#!/bin/bash
#BSUB -J epiDB-Null
#BSUB -o test-%J.out
#BSUB -e test-%J.err
#BSUB -q bigmem
#BSUB -n 1

#Requirements:-
#INSTALL 'seqtk' version 1
#DOWNLOAD 'TE_neoantigens.tsv' and 'TE_chimeric_2297.tsv' from Shah et al. 2023 repository

if [ ! -e "Homo_sapiens_cds_16mers.tab" ] || [ ! -s "Homo_sapiens_cds_16mers.tab" ] ; then
  awk '{if($0~/^>.*/){ln=length(x);for(i=1;i+15<=ln;i++){a[substr(x,i,16)]+=1};x=""}else{x=x $1}}END{for(j in a){print j "\t" a[j]}}' ../Homo_sapiens.GRCh38.cds.all.fa > Homo_sapiens_cds_16mers.tab
fi

awk '(FNR!=1){print $5"\t"$6-1"\t"$7"\t"$1"\t1\t"$11}' TE_chimeric_2297.tsv|sort -Vk1,2 > TE_regions.bed
awk '(NR==FNR){a[$1]=1;next}{OFS="\t";if(a[$4]!=1){$4=$4"x"}; print}' TE_neoantigens.tsv TE_regions.bed > TE_cancerRegions.bed
seqtk subseq -s /data/hemberg/shared_resources/genomes/human/GRCh38.primary_assembly.genome.fa.gz TE_cancerRegions.bed > TE_seqs.fasta

#add TE_id and strand to seqs
awk '(FNR==NR){map[">"$5":"$6"-"$7]=$1"\t"$11;next}{if($0~/^>/){print $1"\t"map[$1]}else{print}}' TE_chimeric_2297.tsv TE_seqs.fasta > TE_seqs_named.fasta
#dustmasker -outfmt fasta -level 10 -in TE_seqs_named.fasta -out TE_seqs_named_masked.fasta
#awk '{if(NR==FNR){a[$1]=1;next}{if($0~/^>.*/){count=0;seqlen=length(seq)-15;for(x=1;x<=seqlen;x++){y=substr(seq,x,16);if(a[y]!=1){print y;count++}};seq=""; print count > "test"}else{seq=seq $1}}}END{count = 0;seqlen=length(seq)-15;for(x=1;x<=seqlen;x++){y=substr(seq,x,16);if(a[y]!=1){print y;count++}};print count > "test"}' ../Homo_sapiens_cds_16mers.tab TE_seqs_named_masked.fasta | grep -v 'a\|t\|g\|c' > TE_nullomers.tsv

awk '{if(NR==FNR){a[$1]=1;next}{if($0~/^>.*/){count=0;seqlen=length(seq)-15;for(x=1;x<=seqlen;x++){y=substr(seq,x,16);if(a[y]!=1){print y;count++}};seq=""; print count > "test"}else{seq=seq $1}}}END{count = 0;seqlen=length(seq)-15;for(x=1;x<=seqlen;x++){y=substr(seq,x,16);if(a[y]!=1){print y;count++}};print count > "test"}' ../Homo_sapiens_cds_16mers.tab TE_seqs_named_masked.fasta > TE_nullomers.all.tsv
#---Keep nullomers with genome-frequency > 300 and overlap < 1 bp
awk -F "\t" '{if(NR==FNR){f[$1]=$2}else{split($1,u,";");c=0;nlist="";olist="";for(x=1;x<length(u);x++){if(f[u[x]]<300){nlist=nlist";"u[x];c++};if(c>6){split(substr(nlist,2),v,";");olist=v[1]";"v[2]";"v[3]";"v[length(v)-2]";"v[length(v)-1]";"v[length(v)]";"}else{olist=substr(nlist,2)}};print olist"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' TE_nullomers_freq.tsv TE_nullomers.all.tsv > TE-nullomersSparse10.tsv

awk '{split($1,u,";");for(i=1;i<length(u);i++)print toupper(u[i])}' ChimerDB-nullomersEdge6.tsv > ChimerDB_nullomers.tsv

#revComp cds_nullomers to search on the 2nd read from RNAseq data
awk 'BEGIN{c["A"]="T";c["T"]="A";c["G"]="C";c["C"]="G";}{y="";for(j=16;j>=1;j--){y=y c[substr($1,j,1)]};print y;}' TE_nullomers.tsv > TE_nullomersRevComp.tsv
