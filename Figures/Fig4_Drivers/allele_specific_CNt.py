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
from matplotlib import font_manager
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines
import seaborn as sns
from matplotlib_venn import venn2
from collections import defaultdict
from collections import Counter
from scipy.stats import norm
from scipy.stats import ttest_ind
from rpy2.robjects import pandas2ri
import rpy2.robjects as ro
from pandas.api.types import CategoricalDtype
from Bio import SeqIO
from PIL import Image
from IPython.display import display
from Bio import Phylo
import ete3
from ete3 import TreeStyle, Tree
from pprint import pprint
import matplotlib.font_manager as fm
import statsmodels.api as sm
import statistics
from scipy import stats
from scipy.stats import f_oneway
import sigProfilerPlotting as sigPlt
from SigProfilerExtractor import subroutines as sub
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from scipy.stats import chisquare
import logging
import itertools
from functools import partial

import inhouse_scripts

chromosomes = [str(i) for i in list(range(1, 23)) + ['X', 'Y', 'MT']]
chrom_cat_type = CategoricalDtype(categories=chromosomes, ordered=True)

seqz_Liver = pd.read_csv("/home/users/wonhee720/Projects/03_CRC/02_vcf/sequenza/52-LRt-O4/52-LRt-O4_segments.txt", header=0, sep="\t") 

seqz_Liver = seqz_Liver[(seqz_Liver['chromosome']=='8') ]


xticks = np.arange(0, 146000000, 30000000)
xtick_labels = ['0','30','60','90',
                '120']


start = seqz_Liver['start.pos']
end = seqz_Liver['end.pos']
start_end = np.ravel(np.array([start,end]), order='F')
mod_start_end = np.reshape(start_end, (-1,2))

A = seqz_Liver.A
A_double = np.array([A, A])
A_double_ravel = np.ravel(A_double, order='F')
mod_A = np.reshape(A_double_ravel, (-1,2))

B = seqz_Liver.B
B_double = np.array([B, B])
B_double_ravel = np.ravel(B_double, order='F')
mod_B = np.reshape(B_double_ravel, (-1,2))

fig, ax = plt.subplots(figsize=(30/25.4,20/25.4))

for pair_x, pair_y in zip(mod_start_end, mod_A):
    if pair_x[1]-pair_x[0] < 1000000 :
        pass
    else:
        ax.plot(pair_x, pair_y, 'm', linewidth=0.5, color='#ff7f0e')
    
for pair_x, pair_y in zip(mod_start_end, mod_B):
    if pair_x[1]-pair_x[0] < 1000000 :
        pass
    else:
        ax.plot(pair_x, pair_y, 'm', linewidth=0.5, color='#1f77b4')

# Plot MYC gene locus (GRCh37 - 8:128747680-128755197)
plt.plot([128751438, 128751438], [0,6.5], linewidth = 0.5, linestyle = '--', color='#e60012')


#plt.grid(True, which='major', linewidth=0.5)
plt.grid(False)

#plt.yticks(ticks=yticks, labels=yticks)z
plt.xticks(ticks=xticks, labels=xtick_labels)
plt.xlim(0,146000000)
plt.yticks(ticks=[0,2,4,6], labels=["0","2","4","6"])


for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(0.5)
    ax.spines[axis].set_color('black')

#ax.tick_params(direction='out', pad=15)
#ax.xaxis.set_tick_params(width=0)
#ax.yaxis.set_tick_params(width=0)
ax.spines['top'].set_color('black')
ax.spines['right'].set_color('black')
ax.spines['left'].set_color('black')
ax.spines['bottom'].set_color('black')

#plt.suptitle("Chr 8",
#             fontproperties='Arial', fontsize = 30, y=1.03)

plt.savefig(f"/home/users/wonhee720/Projects/03_CRC/06_results/drivers_convergent/52-LRt-O4.Chr8.CNt_seqz.svg",format='svg')
plt.show()


# Use germline information to match alleles between samples
MC02_chr8_SNP = scripts.vcf.tidydf("/home/users/wonhee720/Projects/03_CRC/02_vcf/haplotypecaller/hc.snp.merged.readinfo.vep.chr8.vcf.gz").df

dummy_ls = [sampleNameBam(f"01_bam_new/{x}.splitmark.realigned.recal.bam") for x in ["52-B","52-R-FT","52-LRt-O1","52-LRt-O2","52-LRt-O4","52-LRt-O5","52-LRt-O6","52-LRt-O7"]]
dummy_ls = [f"GT:::{x}" for x in dummy_ls]


ls = ["CHROM","POS","REF","ALT"]
MC02_chr8_SNP = MC02_chr8_SNP[ls]

MC02_chr8_SNP["ID"] = '.'
MC02_chr8_SNP["QUAL"] = 0
MC02_chr8_SNP["FILTER"] = "PASS"
MC02_chr8_SNP["AS_FilterStatus"] = "SITE"
MC02_chr8_SNP[dummy_ls] = "."

MC02_chr8_SNP = MC02_chr8_SNP[["CHROM","POS","ID","REF","ALT","QUAL","FILTER","AS_FilterStatus"]+dummy_ls]
MC02_chr8_SNP.rename(columns={"CHROM":"#CHROM"}, inplace=True)

MC02_chr8_SNP.sort_values(by=["#CHROM","POS"], inplace=True)
MC02_chr8_SNP.fillna(".", inplace=True)

display(MC02_chr8_SNP.head())
print(MC02_chr8_SNP.shape)
MC02_chr8_SNP.to_csv("/home/users/wonhee720/Projects/03_CRC/02_vcf/haplotypecaller/hc.snp.merged.readinfo.vep.chr8.modi.tsv", sep="\t",index=False)

# Loading germline chr8 with readinfo annotated
MC02_chr8_SNP = scripts.vcf.tidydf("/home/users/wonhee720/Projects/03_CRC/02_vcf/haplotypecaller/hc.snp.merged.readinfo.vep.chr8.modi.readinfo_tumor.vcf.gz")
MC02_chr8_SNP = MC02_chr8_SNP.df

def homo_hetero_MM_Clip_filter(row, sample):
    sample_VAF = row[f"vaf_all:::{sample}"] 

    # If the SNP is heterozygous, apply filter conditions using reference allele readinfo 
    if (sample_VAF <0.8) & (sample_VAF>0.2):
            sample_threshold = [1,0.05]
        # Filter string when it's a matched normal sample
        #if matched_normal == None:
            #final_sample = (row[f"var_meanMM_all:::{sample}"]<sample_threshold[0])&(row[f"var_ClipFrac_all:::{sample}"] < sample_threshold[1])&(row[f"ref_meanMM_all:::{sample}"]<sample_threshold[0])&(row[f"ref_ClipFrac_all:::{sample}"] < sample_threshold[1])
            final_sample = (row[f"var_meanMM_all:::{sample}"]<sample_threshold[0])&(row[f"var_ClipFrac_all:::{sample}"] < sample_threshold[1])&(row[f"ref_meanMM_all:::{sample}"]<sample_threshold[0])&(row[f"ref_ClipFrac_all:::{sample}"] < sample_threshold[1])
            
            final = final_sample

        # Filter string when it's not a matched normal sample
        #else:
        #    normal_VAF = row[f"vaf_all:::{matched_normal}"]
        #    normal_threshold = [1,0.05] if normal_VAF <0.9 else [2,0.1]
        #    final_sample = (row[f"var_meanMM_all:::{sample}"]<sample_threshold[0])&(row[f"var_ClipFrac_all:::{sample}"] < sample_threshold[1])&(row[f"ref_meanMM_all:::{sample}"]<sample_threshold[0])&(row[f"ref_ClipFrac_all:::{sample}"] < sample_threshold[1])
        #    final_normal = (row[f"var_meanMM_all:::{matched_normal}"]<normal_threshold[0])&(row[f"var_ClipFrac_all:::{matched_normal}"] < normal_threshold[1])&(row[f"ref_meanMM_all:::{matched_normal}"]<normal_threshold[0])&(row[f"ref_ClipFrac_all:::{matched_normal}"] < normal_threshold[1])
        #    final = final_sample & final_normal
    
    # If the SNP is homozygous, DON'T apply filter conditions using reference or variant allele readinfo 
    elif sample_VAF >= 0.8:
            sample_threshold = [2,0.05]
        # Filter string when it's a matched normal sample
        #if matched_normal == None:
            #final_sample = (row[f"var_meanMM_all:::{sample}"]<sample_threshold[0])&(row[f"var_ClipFrac_all:::{sample}"] < sample_threshold[1])
            final_sample = (~np.isnan(row[f"var_meanMM_all:::{sample}"]))&(row[f"var_meanMM_all:::{sample}"]<sample_threshold[0])&(row[f"var_ClipFrac_all:::{sample}"] < sample_threshold[1])
            
            final = final_sample

        # Filter string when it's not a matched normal sample
        #else:
        #    normal_VAF = row[f"vaf_all:::{matched_normal}"]
        #    normal_threshold = [1,0.05] if normal_VAF <0.9 else [2,0.1]
        #    final_sample = (row[f"var_meanMM_all:::{sample}"]<sample_threshold[0])&(row[f"var_ClipFrac_all:::{sample}"] < sample_threshold[1])
        #    final_normal = (row[f"var_meanMM_all:::{matched_normal}"]<normal_threshold[0])&(row[f"var_ClipFrac_all:::{matched_normal}"] < normal_threshold[1])
        #    final = final_sample & final_normal

    elif sample_VAF <= 0.2:
            sample_threshold = [2,0.05]
            final_sample = (~np.isnan(row[f"ref_meanMM_all:::{sample}"]))&(row[f"ref_meanMM_all:::{sample}"]<sample_threshold[0])&(row[f"ref_ClipFrac_all:::{sample}"] < sample_threshold[1])
            final = final_sample
    
    return final


# Divide the total df into per sample dataframes
sample = "52-B"

# Filter string
# Low qualilty reads in target sample
MC02_chr8_SNP = MC02_chr8_SNP[(MC02_chr8_SNP[f"var_read_lowqual:::{sample}"] + MC02_chr8_SNP[f"ref_read_lowqual:::{sample}"] + MC02_chr8_SNP[f"other_read_highqual:::{sample}"] + MC02_chr8_SNP[f"other_read_lowqual:::{sample}"]) 
               < (MC02_chr8_SNP[f"var_read_all:::{sample}"]+MC02_chr8_SNP[f"ref_read_all:::{sample}"] )*0.5]

MC02_chr8_SNP = MC02_chr8_SNP[(MC02_chr8_SNP[f"var_meanMQ_all:::{sample}"]>=60)&(MC02_chr8_SNP[f"var_meanBQ_all:::{sample}"]>=30)]

# Filter calls with excess Mismatch, Clipped portion
MC02_chr8_SNP['MM_Clip_Filter'] = MC02_chr8_SNP.apply(lambda row: homo_hetero_MM_Clip_filter(row, sample), axis=1)
MC02_chr8_SNP = MC02_chr8_SNP[MC02_chr8_SNP['MM_Clip_Filter']]

# Total read count in target sample+
MC02_chr8_SNP = MC02_chr8_SNP[(MC02_chr8_SNP[f"var_read_highqual:::{sample}"] + MC02_chr8_SNP[f"ref_read_highqual:::{sample}"])>10]

MC02_chr8_SNP.to_csv(f"/home/users/wonhee720/Projects/03_CRC/02_vcf/haplotypecaller/{sample}.hc.snp.merged.readinfo.vep.chr8.modi.readinfo_tumor.tsv",
           sep="\t", index=False)




seqz_Primary1 = pd.read_csv("/home/users/wonhee720/Projects/03_CRC/02_vcf/sequenza/52-R-FT/52-R-FT_segments.txt", header=0, sep="\t")


xticks = np.arange(0, 146000000, 30000000)
xtick_labels = ['0','30','60','90',
                '120']

# Remove homozygous SNPs
MC02_chr8_SNP = MC02_chr8_SNP[(MC02_chr8_SNP['vaf_all:::52-B']>0.2) & (MC02_chr8_SNP['vaf_all:::52-B']<0.8)]


# Get heterozygous SNPs whose amplified allele is different between two samples
MC02_chr8_SNP['amp_allele_switched'] = ( (MC02_chr8_SNP['vaf_all:::52-R-FT']>0.5 ) & (MC02_chr8_SNP['vaf_all:::52-LRt-O4']<0.5) ) | ( (MC02_chr8_SNP['vaf_all:::52-R-FT']<0.5 ) & (MC02_chr8_SNP['vaf_all:::52-LRt-O4']>0.5) ) 

start = seqz_Primary1['start.pos']
end = seqz_Primary1['end.pos']
start_end = np.ravel(np.array([start,end]), order='F')
mod_start_end = np.reshape(start_end, (-1,2))

A = seqz_Primary1.A
A_double = np.array([A, A])
A_double_ravel = np.ravel(A_double, order='F')
mod_A = np.reshape(A_double_ravel, (-1,2))

B = seqz_Primary1.B
B_double = np.array([B, B])
B_double_ravel = np.ravel(B_double, order='F')
mod_B = np.reshape(B_double_ravel, (-1,2))

fig, ax = plt.subplots(figsize=(30/25.4,20/25.4))

for pair_x, pair_y in zip(mod_start_end, mod_A):
    if pair_x[1]-pair_x[0] < 1000000 :
        pass
    else:
        start_pos = pair_x[0]
        end_pos = pair_x[1]
        _df = MC02_chr8_SNP[(MC02_chr8_SNP['POS']>start_pos) & (MC02_chr8_SNP['POS']<end_pos)]
        
        if (_df.shape[0]*0.75)<_df['amp_allele_switched'].sum(axis=0):
            ax.plot(pair_x, pair_y, 'm', linewidth=0.5, color='#1f77b4')
        else:
            ax.plot(pair_x, pair_y, 'm', linewidth=0.5, color='#ff7f0e')
    
for pair_x, pair_y in zip(mod_start_end, mod_B):
    if pair_x[1]-pair_x[0] < 1000000 :
        pass
    else:
        start_pos = pair_x[0]
        end_pos = pair_x[1]
        _df = MC02_chr8_SNP[(MC02_chr8_SNP['POS']>start_pos) & (MC02_chr8_SNP['POS']<end_pos)]
        
        if (_df.shape[0]*0.75)<_df['amp_allele_switched'].sum(axis=0):
            ax.plot(pair_x, pair_y, 'm', linewidth=0.5, color='#ff7f0e')
        else:
            ax.plot(pair_x, pair_y, 'm', linewidth=0.5, color='#1f77b4')

# Plot MYC gene locus (GRCh37 - 8:128747680-128755197)
plt.plot([128751438, 128751438], [0,8.5], linewidth = 0.5, linestyle = '--', color='#e60012')


#plt.grid(True, which='major', linewidth=0.5)
plt.grid(False)

#plt.yticks(ticks=yticks, labels=yticks)z
plt.xticks(ticks=xticks, labels=xtick_labels)
plt.xlim(0,146000000)
plt.ylim((0,10))
plt.yticks(ticks=[0,2,4,6,8], labels=["0","2","4","6","8"])


for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(0.5)
    ax.spines[axis].set_color('black')

#ax.tick_params(direction='out', pad=15)
#ax.xaxis.set_tick_params(width=0)
#ax.yaxis.set_tick_params(width=0)
ax.spines['top'].set_color('black')
ax.spines['right'].set_color('black')
ax.spines['left'].set_color('black')
ax.spines['bottom'].set_color('black')

#plt.suptitle("Chr 8",
#             fontproperties='Arial', fontsize = 30, y=1.03)

plt.savefig(f"/home/users/wonhee720/Projects/03_CRC/06_results/drivers_convergent/52-R-FT.Chr8.CNt_seqz.svg",format='svg')
plt.show()
