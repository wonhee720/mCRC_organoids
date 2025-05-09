#!/bin/bash
# $1 : sample id
# $2 : chromosome
# $3 : start pos
# $4 : end pos
# $5 : cluster number

/home/users/wonhee720/anaconda3/envs/genome_new/bin/Rscript \
/mCRC/Figures/Fig4_Drivers/yilong_plot_draw.R \
"/home/users/wonhee720/Projects/03_CRC/02_vcf/mosdepth/$1/$1.cov_corrected.abs.tsv.absCN.gen_fi" \
/home/users/wonhee720/Projects/03_CRC/02_vcf/gremlin_rescue/apply/$1/$1.gremlin.somatic.svs.concat.bptinfo.anv.filtered_1st.tsv.clustered.gap_seg.100kbAbsCN.complex_class.branch.cluster$5 \
$2 \
$3 \
$4 \
$5