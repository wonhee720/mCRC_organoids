#!/bin/bash

set -eu

cd ~/Projects/03_CRC/02_vcf/union

tsv=$(find ./ -type f -name "*$1*.filtered_1st.tsv")

for i in ${tsv}
do
name=$(echo ${i}| cut -d'/' -f2)
echo $name
seq="/home/users/wonhee720/Projects/03_CRC/02_vcf/sequenza/"$name"/"$name"_segments.txt"
~/anaconda3/envs/genome_new/bin/python "/mCRC/inhouse_scripts/SNV_CleanSeqz_CN_annot.py" \
${i} ${seq}
done
