import os
import sys
import re
import copy
import yaml
import array
import warnings
import pysam 
import cyvcf2
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
from Bio import Phylo
from scipy import stats
from scipy.stats import f_oneway
from scipy.stats import chisquare
from scipy.stats import ttest_ind
import io
import ete3
import math
import logging
import matplotlib.font_manager as fm
import matplotlib.ticker as ticker
import itertools
from matplotlib.ticker import MaxNLocator
from scipy.stats import chi2_contingency

import inhouse_scripts

chromosomes = [str(i) for i in list(range(1, 20)) + ['X', 'Y']]
chrom_cat_type = CategoricalDtype(categories=chromosomes, ordered=True)



df1 = pd.read_csv(path1, sep="\t")
df2 = pd.read_csv(path2, sep="\t")


# Get a function to annotate proportion of reads that have mates mapped to other chromosomes
# Use this column to filter out SVs

def get_mate_chrs(samfile, row):
    CHR1 = row["#CHR1"]
    POS1 = row["POS1"]
    CHR2 = row["CHR2"]
    POS2 = row["POS2"]

    count_same1 = 0 # Number of mate reads at the same chromosome
    count_diff1 = 0 # Number of mate reads at the different chromosome
    count_same2 = 0
    count_diff2 = 0
    
    # Get 1kb area of the breakpoint
    iterr1 = samfile.fetch(CHR1, POS1-base_range, POS1+base_range)
    iterr2 = samfile.fetch(CHR2, POS2-base_range, POS2+base_range)
    
    for x in iterr1:
        mate_chr = str(x).split("\t")[6]
        if mate_chr == '22':
            mate_chr = 'X' 
        elif mate_chr == '23':
            mate_chr = 'Y'
        elif int(mate_chr) < 22:
            mate_chr = str(int(mate_chr)+1)
        else:
            pass
        
        if mate_chr == str(CHR1):
            count_same1 += 1
        else:
            count_diff1 += 1
            
    prop1 = count_diff1/(count_same1+count_diff1) # The proportion of mate reads at the different chromosome
    
    for x in iterr2:
        mate_chr = str(x).split("\t")[6]
        if mate_chr == '22':
            mate_chr = 'X' 
        elif mate_chr == '23':
            mate_chr = 'Y'
        elif int(mate_chr) < 22:
            mate_chr = str(int(mate_chr)+1)
        else:
            pass
        
        if mate_chr == str(CHR2):
            count_same += 1
        else:
            count_diff += 1
            
    prop2 = count_diff2/(count_same2+count_diff2)
        
    return count_same1, count_diff1, prop1, count_same2, count_diff2, prop2


# Annotate
samfile = pysam.AlignmentFile(bam_path, "rb")
for __, row in df.iterrows():
    get_mate_chrs(samfile, row)


# Filtering tier 1 (df1)
df1 = df1[ (df1["normal_var_read"]<1) & (df1['normal_var_split_read']<1) & (df1['normal_same_clip']<1) 
          & (df1["normal_panel"]<1) & (df1['tumor_var_mapq_bpt1'] == 60) & (df1['score'] > 0.95) 
          & ((df1['normal_other_var_cluster_bpt1']<5) & (df1['normal_other_var_cluster_bpt2']<5))] # Filter if any var read exist in normal sample, PoN > 0, poor MQ area, with gremlin score, other variant clusters >=5 in normal in one of breakpoint


# Split the dataframe according to SV type
df1_inv = df1[df1["SVTYPE"] == 'INV'] # INV
df1_tra = df1[df1["SVTYPE"] == 'TRA'] # TRA
df1_rest = df1[(df1["SVTYPE"] == 'DUP') | (df1["SVTYPE"] == 'DEL') ] # DEL, DUP

# Filtering for INV
df1_inv = df1_inv[(df1_inv["tumor_var_read"]>20) | (df1_inv["tumor_var_split_read"]>3)] # For INV, apply more strict conditions -> var_read > 20 or var_split read >3
df1_inv = df1_inv[((df1_inv["tumor_var_split_read"]>10) & (df1_inv['mate_chr_diff_prop_bp1']<0.3) & (df1_inv['mate_chr_diff_prop_bp2']<0.3)) | 
    (df1_inv["tumor_var_split_read"]<=10) & (df1_inv['mate_chr_diff_prop_bp1']<0.1) & (df1_inv['mate_chr_diff_prop_bp2']<0.1)] # if split read >10, apply more generous mate_chr_diff_prop threshold, if split read <=10, apply more strict mate_chr_diff_prop threshold

# Filtering for DEL, DUP
df1_rest = df1_rest[(df1_rest["tumor_var_read"]>10) | (df1_rest["tumor_var_split_read"]>0)] # For DUP, DEL, apply more generous conditions -> var_read > 10 or var_split read > 0
df1_rest = df1_rest[((df1_rest["tumor_var_split_read"]>3) & (df1_rest['mate_chr_diff_prop_bp1']<0.3) & (df1_rest['mate_chr_diff_prop_bp2']<0.3)) | 
    (df1_rest["tumor_var_split_read"]<=3) & (df1_rest['mate_chr_diff_prop_bp1']<0.1) & (df1_rest['mate_chr_diff_prop_bp2']<0.1)] # if split read >3, apply more generous mate_chr_diff_prop threshold, if split read <=3, apply more strict mate_chr_diff_prop threshold


# Filtering for TRA
df1_tra = df1_tra[(df1_tra["tumor_var_read"]>10) | (df1_tra["tumor_var_split_read"]>0)]


# Concat the split dataframes into one once again
df1 = pd.concat([df1_inv,_df_rest,_df_tra])



# Filtering tier 2 (df2)
df2 = df2[ (df2["normal_var_read"]<1) & (df2['normal_var_split_read']<1) & (df2['normal_same_clip']<1) & 
      (df2["normal_panel"]<1) & (df2['tumor_var_mapq_bpt1'] == 60) & (df2['score'] > 0.5) &
      ((df2['normal_other_var_cluster_bpt1']<5) & (df2['normal_other_var_cluster_bpt2']<5))] # Filter if any var read exist in normal sample, PoN > 0, poor MQ area, with gremlin score, other variant clusters >=5 in normal in one of breakpoint


# Additional filtering for tier 2 compared to tier 1
df2 = df2[(df2['tumor_var_split_read']>0) & (df2['tumor_var_read']>3)]


# Split the dataframe according to SV type
df2_inv = df2[df2["SVTYPE"] == 'INV'] # INV
df2_tra = df2[df2["SVTYPE"] == 'TRA'] # TRA
df2_rest = df2[(df2["SVTYPE"] == 'DUP') | (df2["SVTYPE"] == 'DEL') ] # DEL, DUP


# Filtering for INV - More strict filtering than tier 1
df2_inv['distance'] = abs(df2_inv['POS1'] - df2_inv['POS2'])
df2_inv = df2_inv[~( (df2_inv['distance']<10000) & (df2_inv['tumor_var_split_read']<5) ) ]
df2_inv = df2_inv[(df2_inv["tumor_var_read"]>20) | (df2_inv["tumor_var_split_read"]>3)] # For INV, apply more strict conditions -> var_read > 20 or var_split read >3
df2_inv = df2_inv[((df2_inv["tumor_var_split_read"]>10) & (df2_inv['mate_chr_diff_prop_bp1']<0.3) & (df2_inv['mate_chr_diff_prop_bp2']<0.3)) | 
(df2_inv["tumor_var_split_read"]<=10) & (df2_inv['mate_chr_diff_prop_bp1']<0.1) & (df2_inv['mate_chr_diff_prop_bp2']<0.1)] # if split read >10, apply more generous mate_chr_diff_prop threshold, if split read <=10, apply more strict mate_chr_diff_prop threshold


# Filtering for DEL, DUP
df2_rest = df2_rest[(df2_rest["tumor_var_read"]>10) | (df2_rest["tumor_var_split_read"]>0)] # For DUP, DEL, apply more generous conditions -> var_read > 10 or var_split read > 0
df2_rest = df2_rest[((df2_rest["tumor_var_split_read"]>3) & (df2_rest['mate_chr_diff_prop_bp1']<0.3) & (df2_rest['mate_chr_diff_prop_bp2']<0.3)) | 
(df2_rest["tumor_var_split_read"]<=3) & (df2_rest['mate_chr_diff_prop_bp1']<0.1) & (df2_rest['mate_chr_diff_prop_bp2']<0.1)] # if split read >3, apply more generous mate_chr_diff_prop threshold, if split read <=3, apply more strict mate_chr_diff_prop threshold


# Filtering for TRA
df2_tra = df2_tra[(df2_tra["tumor_var_read"]>10) | (df2_tra["tumor_var_split_read"]>0)]


# Concat the split dataframes into one once again
df2_inv = df2_inv.drop('distance', axis=1, inplace=True)
df2 = pd.concat([df2_inv,df2_rest,df2_tra])


concat_df = pd.concat([df1,df2])
concat_df.to_csv(path, sep="\t", index=False)