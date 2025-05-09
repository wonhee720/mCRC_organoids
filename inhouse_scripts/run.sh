#!/bin/bash
set -eu


prog=/home/users/pjh/scripts/annotation/short_variants/readinfoAnnot.wrapper.v2.sh


nbam=/home/users/kimin/projects/03_Mucosal_Melanoma/01_bam/MM001-Blood.splitmark.realigned.recal.bam
tbam=/home/users/kimin/projects/03_Mucosal_Melanoma/01_bam/MM001-Tumor.splitmark.realigned.recal.bam
nid=normal
tid=tumor

in_vcf=merged.snv.vcf.gz
out_vcf=merged.snv.readinfo.vcf.gz

fasta=/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta

$prog \
	-i $in_vcf \
	-b ${nbam},${tbam} \
	-s ${nid},${tid} \
	--sched local \
	-f $fasta \
	-o $out_vcf \
	-O z

## VEP ##

vep_prog=/home/users/pjh/scripts/annotation/short_variants/vep_annotation_wrapper.v6.sh

in_vcf2=$out_vcf
out_vcf2=merged.snv.readinfo.vep.vcf.gz

$vep_prog \
	-o $out_vcf2 \
	-O z \
	-v 19 \
	-p 1 \
	--sched local \
	$in_vcf2


## pon ##
pon_snv=/home/users/pjh/References/PON19/merged_snv/pon_merged.hg19.snuN30.pcawgN36.bgiN24.pcnslN11.snv.bcf.gz
pon_indel=/home/users/pjh/References/PON19/merged_indel/pon_merged.hg19.snuN30.pcawgN36.bgiN24.pcnslN11.indel.bcf.gz

in_vcf3=$out_vcf2
out_vcf3=merged.snv.readinfo.vep.pon.vcf.gz

bcftools annotate $in_vcf3 -a $pon_snv -c INFO -O z -o $out_vcf3
bcftools index $out_vcf3


## dbsnp ##
dbsnp_db=/home/users/pjh/References/dbSNP37/modified_files_201214/dbSNP_b154_GRCh37.p13.bcf.gz

in_vcf4=$out_vcf3
out_vcf4=merged.snv.readinfo.vep.pon.dbsnp.vcf.gz

bcftools annotate $in_vcf4 -a $dbsnp_db -c INFO/dbSNP_MAX_FREQ,INFO/dbSNP_MAX_FREQ_POP -O z -o $out_vcf4
bcftools index $out_vcf4
