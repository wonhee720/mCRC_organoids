# 20210409 created by kimin
# reheadering is important

input_bam=$1
output_bam=$2

samtools view -H $input_bam|grep -v Mus | samtools reheader - $input_bam > $output_bam
samtools index $output_bam