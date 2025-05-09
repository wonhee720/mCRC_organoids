# Script for preparing dataframes for R circlize
# modified to receive gc corrected mosdepth for cnv - 22.02.13

import pandas as pd
import numpy as np
import os, sys
import pysam

sys.path.append('/home/users/wonhee720/Projects/03_CRC/')

import scripts

fasta = pysam.FastaFile("/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta")
chromosomes = [str(i) for i in list(range(1, 23)) + ['X', 'Y', 'MT']]


sample = sys.argv[1]
path_snv = sys.argv[2]
path_indel = sys.argv[3]
path_cnv = sys.argv[4]
path_sv = sys.argv[5]

path_snv_out = sys.argv[6]
path_indel_out = sys.argv[7]
path_cnv_out = sys.argv[8]
path_sv_out = sys.argv[9]


def sampleNameBam(bamFile):
    """get @RG SM: information as sample name from BAM header"""
    bam = pysam.AlignmentFile(bamFile)
    name = bam.header['RG'][0]['SM']
    return name


# Prepare SNV
snv = pd.read_csv(path_snv, sep="\t", dtype={"#CHROM": str})
tmp_dfs = []

for chrom in chromosomes:
    
    df_tmp = snv[snv["#CHROM"] == str(chrom)]
    df_tmp['DIFF'] = df_tmp['POS'].diff()
    df_tmp['logDIFF'] = np.log10(df_tmp['DIFF'])
    tmp_dfs.append(df_tmp)
    
snv = pd.concat(tmp_dfs)
# Need to make context_3 column
snv[['context_3', 'substitution']] = snv.query("VAR_TYPE == 'snp'").apply(lambda x: scripts.sig.context_3(fasta, x['#CHROM'], x['POS'], x['ALT']), axis=1, result_type='expand')

snv = snv[["#CHROM", "POS","logDIFF", "context_3"]]
snv["#CHROM"] = "chr" + snv["#CHROM"].astype(str)
snv.rename(columns = {"#CHROM":"chr"}, inplace=True)
snv = snv[["chr","POS","logDIFF","context_3"]]
snv = snv.fillna(0)
snv["type"] = snv["context_3"].str.slice(start=2, stop=5)

snv.to_csv(path_snv_out, sep="\t", index=False)


# Prepare Indels
indel = pd.read_csv(path_indel, sep="\t", encoding='ISO-8859-1')
tumor_original = sampleNameBam(f"/home/users/wonhee720/Projects/03_CRC/01_bam_new/{sample}.splitmark.realigned.recal.bam")

indel = indel[["#CHROM", "POS", "REF", "ALT", f"vaf_all:::{tumor_original}"]] 
indel["#CHROM"] = "chr" + indel["#CHROM"].astype(str)
indel.rename(columns={"#CHROM":"chr", f"vaf_all:::{tumor_original}":"vaf"}, inplace=True)
indel["type"] = indel["REF"] > indel["ALT"]

indel["type"].replace({True:"del",False:"ins"}, inplace=True)

indel.to_csv(path_indel_out, sep="\t", index=False)


# Prepare CNV
###### CNV using sequenza ######
#cnv = pd.read_csv(path_cnv, sep="\t")
#cnv = cnv[["chromosome", "start.pos","end.pos","CNt","A","B"]]
#cnv.rename(columns={"chromosome":"chr", "start.pos":"start","end.pos":"end"}, inplace=True)
#cnv["chr"] = "chr" + cnv["chr"].astype(str)
#cnv['diff'] = abs(cnv['start'] - cnv['end'])
#cnv = cnv[cnv['diff']>100000] # Remove segments smaller than 100kbp
#cnv.to_csv(path_cnv_out,sep="\t", index=False)


###### CNV using gc corrected (abs) mosdepth ######
cnv = pd.read_csv(path_cnv, sep="\t", header=0, dtype={'chr':str,'pos':int,'pos_2':int,'GC':float,'average_cov':float,'cov_corrected':float}, encoding="utf-8")
cnv = cnv[["chr", "pos","cov_corrected_abs"]]
cnv.rename(columns={"pos":"start"}, inplace=True)
cnv["chr"] = "chr" + cnv["chr"].astype(str)
cnv.to_csv(path_cnv_out, sep="\t", index=False, encoding="utf-8")


# Prepare SV
sv = pd.read_csv(path_sv, sep="\t")
sv.rename(columns={"#CHR1":"CHR1"}, inplace=True)
sv["CHR1"] = "chr" + sv["CHR1"].astype(str)
sv["CHR2"] = "chr" + sv["CHR2"].astype(str)
sv=sv[["CHR1","POS1","CHR2","POS2","CT","SVTYPE"]]
sv.to_csv(path_sv_out,sep="\t", index=False)