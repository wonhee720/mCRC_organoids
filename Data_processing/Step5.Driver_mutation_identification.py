import os
import sys
import re
import copy
import yaml
import array
import warnings
import pysam 
import cyvcf2
import copy
# import igv
import allel
import time
import tempfile
import subprocess
import shlex
import pyranges as pr
import pickle as pkl
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2
from collections import defaultdict
from collections import Counter
from scipy.stats import norm
from rpy2.robjects import pandas2ri
import rpy2.robjects as ro
from pandas.api.types import CategoricalDtype
from Bio import SeqIO
from PIL import Image
from IPython.display import display

import inhouse_scripts

chromosomes = [str(i) for i in list(range(1, 23)) + ['X', 'Y', 'MT']]
chrom_cat_type = CategoricalDtype(categories=chromosomes, ordered=True)

tcga = pd.read_csv("/mCRC/Data_preprocessing/cancer_gene_census_Cosmic.csv")


snv_filtered_1st_tsv_path_dict = {k : [f"/home/users/wonhee720/Projects/03_CRC/02_vcf/union/{x}/{x}.snps.filtered_1st.tsv-seqzcn-ccf" for x in sample_id_dict[k]] for k in patient_id_list }
snv_filtered_1st_vcf_path_dict = {k : [f"/home/users/wonhee720/Projects/03_CRC/02_vcf/union/{x}/{x}.snps.filtered_1st.seqzcn-ccf.vcf.gz" for x in sample_id_dict[k]] for k in patient_id_list }

indel_filtered_1st_tsv_path_dict = {k : [f"/home/users/wonhee720/Projects/03_CRC/02_vcf/union/{x}/{x}.indels.filtered_1st.tsv-seqzcn-ccf" for x in sample_id_dict[k]] for k in patient_id_list }
indel_filtered_1st_vcf_path_dict = {k : [f"/home/users/wonhee720/Projects/03_CRC/02_vcf/union/{x}/{x}.indels.filtered_1st.seqzcn-ccf.vcf.gz" for x in sample_id_dict[k]] for k in patient_id_list }

mnv_filtered_1st_tsv_path_dict = {k : [f"/home/users/wonhee720/Projects/03_CRC/02_vcf/union/{x}/{x}.mnps.filtered_1st.tsv-seqzcn-ccf" for x in sample_id_dict[k]] for k in patient_id_list }
mnv_filtered_1st_vcf_path_dict = {k : [f"/home/users/wonhee720/Projects/03_CRC/02_vcf/union/{x}/{x}.mnps.filtered_1st.seqzcn-ccf.vcf.gz" for x in sample_id_dict[k]] for k in patient_id_list }

sequenza_path_dict = {k : [f"/home/users/wonhee720/Projects/03_CRC/02_vcf/sequenza/{x}/{x}_segments.txt" for x in sample_id_dict[k]] for k in patient_id_list }

sv_path_dict = {k : [f"/home/users/wonhee720/Projects/03_CRC/02_vcf/gremlin_rescue/apply/{x}.gremlin.somatic.svs.concat.bptinfo.anv.filtered_1st.tsv" for x in sample_id_dict[k]] for k in patient_id_list }


# SNV
dfs_snv_filtered_1st = {}

for k,v in snv_filtered_1st_vcf_path_dict.items():
    print("patient " + str(k) + " !!!")
    dfs_snv_filtered_1st[k] = {}
    for sample, path in zip(sample_id_dict[k], v):
        dfs_snv_filtered_1st[k][sample] = scripts.vcf.tidydf(path)
        print(path + " loaded !!!")

# Set dictionary for snv, indel, mnv driver candidate concatenation
driver_concat_dfs = defaultdict(dict)

for patient in dfs_snv_filtered_1st.keys():
    driver_concat_dfs[patient] = defaultdict(list)

# Searching among SNVs
for patient in dfs_snv_filtered_1st.keys():
    for k, v in dfs_snv_filtered_1st[patient].items():
        
        _df = v.df
        _df["CSQ_SYMBOL_canonical"] = _df["CSQ_SYMBOL_canonical"].str.split(",")
        _df = _df.explode("CSQ_SYMBOL_canonical")
        _df_merged =_df.merge(tcga,how="left", left_on="CSQ_SYMBOL_canonical", right_on="Gene Symbol", suffixes=('_sample','_tcga'))
        _df_drivers = _df_merged[~_df_merged["Gene Symbol"].isna()]
        _df_drivers = _df_drivers[~_df_drivers["CSQ_IMPACT_canonical"].str.contains("LOW")]
        _df_drivers = _df_drivers[_df_drivers["CSQ_IMPACT_canonical"].str.contains("HIGH")|_df_drivers["CSQ_IMPACT_canonical"].str.contains("MODERATE")|
                                  _df_drivers["CSQ_Consequence_canonical"].str.contains("coding_sequence_variant|mature_miRNA_variant|NMD_transcript_variant|TFBS_ablation|TFBS_amplification|TF_binding_site_variant|regulatory_region_ablation|regulatory_region_amplification|regulatory_region_variant", regex=True)]
    
        number = _df_drivers.shape[0]
        driver_concat_dfs[patient][k].append(_df_drivers)
        # print(f"{tumor} has {number} possible driver mutations !!!")
        # display(_df_drivers[["CHROM","POS","REF","ALT",f"vaf:::{tumor_original}","Gene Symbol","Role in Cancer","CSQ_Consequence_canonical","CSQ_IMPACT_canonical","Tumour Types(Somatic)","COSMIC_ID","COSMIC_occurrence","COSMIC_total_occurrence"]])


# Indel

dfs_indel_filtered_1st = {}

for k,v in indel_filtered_1st_vcf_path_dict.items():
    print("patient " + str(k) + " !!!")
    dfs_indel_filtered_1st[k] = {}
    for sample, path in zip(sample_id_dict[k], v):
        dfs_indel_filtered_1st[k][sample] = scripts.vcf.tidydf(path)
        print(path + " loaded !!!")

# Searching among Indels
for patient in dfs_indel_filtered_1st.keys():
    for k, v in dfs_indel_filtered_1st[patient].items():
        _df = v.df
        _df["CSQ_SYMBOL_canonical"] = _df["CSQ_SYMBOL_canonical"].str.split(",")
        _df = _df.explode("CSQ_SYMBOL_canonical")
        _df_merged =_df.merge(tcga,how="left", left_on="CSQ_SYMBOL_canonical", right_on="Gene Symbol", suffixes=('_sample','_tcga'))
        _df_drivers = _df_merged[~_df_merged["Gene Symbol"].isna()]
        _df_drivers = _df_drivers[~_df_drivers["CSQ_IMPACT_canonical"].str.contains("LOW")]
        _df_drivers = _df_drivers[_df_drivers["CSQ_IMPACT_canonical"].str.contains("HIGH")|_df_drivers["CSQ_IMPACT_canonical"].str.contains("MODERATE")|
                                  _df_drivers["CSQ_Consequence_canonical"].str.contains("coding_sequence_variant|mature_miRNA_variant|NMD_transcript_variant|TFBS_ablation|TFBS_amplification|TF_binding_site_variant|regulatory_region_ablation|regulatory_region_amplification|regulatory_region_variant", regex=True)]
        
        number = _df_drivers.shape[0]
        driver_concat_dfs[patient][k].append(_df_drivers)
        # print(f"{tumor} has {number} possible driver mutations !!!")
        # display(_df_drivers[["CHROM","POS","REF","ALT",f"vaf:::{tumor_original}","vaf:::TUMOR","Gene Symbol","Role in Cancer","CSQ_Consequence_canonical","CSQ_IMPACT_canonical","Tumour Types(Somatic)","COSMIC_ID","COSMIC_occurrence","COSMIC_total_occurrence"]])


# MNV
dfs_mnv_filtered_1st = {}

for k,v in mnv_filtered_1st_vcf_path_dict.items():
    print("patient " + str(k) + " !!!")
    dfs_mnv_filtered_1st[k] = {}
    for sample, path in zip(sample_id_dict[k], v):
        dfs_mnv_filtered_1st[k][sample] = scripts.vcf.tidydf(path)
        print(path + " loaded !!!")


# Searching among MNVs
for patient in dfs_mnv_filtered_1st.keys():
    for k, v in dfs_mnv_filtered_1st[patient].items():
    
        _df = v.df
        _df["CSQ_SYMBOL_canonical"] = _df["CSQ_SYMBOL_canonical"].str.split(",")
        _df = _df.explode("CSQ_SYMBOL_canonical")
        _df_merged =_df.merge(tcga,how="left", left_on="CSQ_SYMBOL_canonical", right_on="Gene Symbol", suffixes=('_sample','_tcga'))
        _df_drivers = _df_merged[~_df_merged["Gene Symbol"].isna()]
        _df_drivers = _df_drivers[~_df_drivers["CSQ_IMPACT_canonical"].str.contains("LOW")]
        _df_drivers = _df_drivers[_df_drivers["CSQ_IMPACT_canonical"].str.contains("HIGH")|_df_drivers["CSQ_IMPACT_canonical"].str.contains("MODERATE")|
                                  _df_drivers["CSQ_Consequence_canonical"].str.contains("coding_sequence_variant|mature_miRNA_variant|NMD_transcript_variant|TFBS_ablation|TFBS_amplification|TF_binding_site_variant|regulatory_region_ablation|regulatory_region_amplification|regulatory_region_variant", regex=True)]
        
        number = _df_drivers.shape[0]
        driver_concat_dfs[patient][k].append(_df_drivers)
        # print(f"{tumor} has {number} possible driver mutations !!!")
        # display(_df_drivers[["CHROM","POS","REF","ALT",f"vaf:::{tumor_original}","vaf:::TUMOR","Gene Symbol","Role in Cancer","CSQ_Consequence_canonical","CSQ_IMPACT_canonical","Tumour Types(Somatic)","COSMIC_ID","COSMIC_occurrence","COSMIC_total_occurrence"]])

# Concat SNV, Indel, MNV driver candidates
driver_final_dfs = defaultdict(dict)

for patient in dfs_mnv_filtered_1st.keys():
    for k, v in driver_concat_dfs[patient].items():
        sum_df = pd.concat(v)
        driver_final_dfs[patient][k] = sum_df


for patient in dfs_mnv_filtered_1st.keys():
    for (k, v), seqz in zip(driver_final_dfs[patient].items(), sequenza_path_dict[patient]):
        tumor = sampleNameBam(f"01_bam_new/{k}.splitmark.realigned.recal.bam")
        
        seqz_df = pd.read_csv(seqz, sep="\t")
        v = v.rename(columns = {"CHROM":"#CHROM"})
        v = v.reset_index(drop=True)
        rows = v.shape[0]
        v["CNt"] = np.nan
        v["A"] = np.nan
        v["B"] = np.nan
    
        for i in range(rows):
            chrom = v.loc[i,"#CHROM"]
            pos = v.loc[i,"POS"]
            cnv = seqz_df[((seqz_df["chromosome"] == chrom) & (seqz_df["start.pos"] < pos) & (seqz_df["end.pos"] > pos))][["CNt","A","B"]]
    
            try:
                CNt = cnv.iloc[0,0]
                A = cnv.iloc[0,1]
                B = cnv.iloc[0,2]
                v.loc[i,"CNt"] = CNt
                v.loc[i,"A"] = A
                v.loc[i,"B"] = B
        
            except IndexError:
                pass
        
        print(k)
        driver_final_dfs[patient][k] = v
        
        display(v[["#CHROM","POS","REF","ALT","CNt","A","B",f"vaf_all:::{tumor}","Gene Symbol","Role in Cancer","CSQ_Consequence_canonical","CSQ_CLIN_SIG_canonical","CSQ_IMPACT_canonical","COSMIC_occurrence","COSMIC_total_occurrence"]])
        v[["#CHROM","POS","REF","ALT","CNt","A","B",f"vaf_all:::{tumor}","Gene Symbol","Role in Cancer","CSQ_Consequence_canonical","CSQ_CLIN_SIG_canonical","CSQ_IMPACT_canonical","COSMIC_occurrence","COSMIC_total_occurrence"]].to_csv(f"06_results/drivers/{k}.drivers.tsv", sep="\t", index=False)





















