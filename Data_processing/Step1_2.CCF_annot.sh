#!/bin/bash

set -eu

cd ~/Projects/03_CRC/02_vcf/union

tsv=$(find ./ -type f -name "*$1*.filtered_1st.tsv-seqzcn")

for i in ${tsv}
do
name=$(echo ${i}|cut -d'/' -f2)
echo $name
seq="/home/users/wonhee720/Projects/03_CRC/02_vcf/sequenza/"$name"/"$name"_confints_CP.txt"
~/anaconda3/envs/genome_new/bin/python /home/users/wonhee720/Projects/03_CRC/scripts_lwh/SNV_mutCN_CCF_annot.CRC.py \
${i} ${seq} ${name}
done