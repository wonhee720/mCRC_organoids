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



mutsig_LINE1_df = pd.read_csv("/home/users/wonhee720/Projects/03_CRC/06_results/mutcounts.sig_snv_indel_per_sample.LINE1.from_branch_concat.240130_fixed.tsv",
                             sep="\t", index_col=0, header=0,
                             dtype={'Patient':str})

mutsig_LINE1_df.reset_index(names="sampleID", inplace=True)


# Some modifications to concat the dataframe
PCAWG_SBS_CRCmodi_L1['Sample'] = "Cancer"
PCAWG_SBS_CRCmodi_L1['Chemotherapy duration'] = 0
PCAWG_SBS_CRCmodi_L1["sampleID"] = PCAWG_SBS_CRCmodi_L1['Sample Names']
PCAWG_SBS_CRCmodi_L1['Patient'] = None

NormalColon_L1_df['Sample'] = "Normal"
NormalColon_L1_df['Chemotherapy duration'] = 0
NormalColon_L1_df['Patient'] = None

mutsig_LINE1_df['Sample'] = 'Cancer organoid'


important_cols = ["sampleID","Sample","Chemotherapy duration","L1 rate","Total L1",'Patient']
Concat_cohorts_L1rate = pd.concat([NormalColon_L1_df[important_cols], PCAWG_SBS_CRCmodi_L1[important_cols], mutsig_LINE1_df[important_cols]])


Concat_cohorts_L1rate["Patient"] = Concat_cohorts_L1rate["Patient"].apply(lambda x : pt_id_exchanger[x] if x!=None else x)

Concat_cohorts_L1rate['Sample2'] = Concat_cohorts_L1rate['Sample'] + "_" + Concat_cohorts_L1rate['Chemotherapy duration'].astype(str)
Concat_cohorts_L1rate['Sample3'] = Concat_cohorts_L1rate['Sample'] + "_" + Concat_cohorts_L1rate['Patient'].astype(str)

pt_color_dict_modi2 = {"Normal_None": "#c8c8c8",
'Cancer_None': "#646464" ,
'Cancer organoid_MC05': '#1F77B4',
 'Cancer organoid_MC06': '#FF7F0E',
 'Cancer organoid_MC01': '#8C564B',
 'Cancer organoid_MC02': '#9467BD',
 'Cancer organoid_MC04': '#2CA02C',
 'Cancer organoid_MC03': '#D62728'}

 # Remove FT samples
Concat_cohorts_L1rate = Concat_cohorts_L1rate[~Concat_cohorts_L1rate['sampleID'].str.contains("-F")]


# Scatterplot of L1 rate per EPM (endogenous point mutations) per sample (grouped by chemotherapy durations)
fig, ax = plt.subplots(figsize=(38/25.4 , 30/25.4))
#fig, ax = plt.subplots(figsize=(15 , 15))

sns.boxplot(ax=ax, data=Concat_cohorts_L1rate, x='Sample2', y='Total L1', color='white', linewidth=0.5, linecolor='black', fliersize=0, legend=False)
sns.stripplot(ax=ax, data=Concat_cohorts_L1rate, x='Sample2', y='Total L1', hue='Sample3', palette=pt_color_dict_modi2,
              dodge=False, linewidth=0.15, size=2, zorder=6, legend=False, edgecolor='black')

plt.setp(ax.patches, linewidth=0.5, edgecolor='black')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)
ax.spines['left'].set_linewidth(0.5)
ax.spines['bottom'].set_linewidth(0.5)
ax.spines['left'].set_color('black')
ax.spines['bottom'].set_color('black')

ax.set_xticklabels(['NL-O', 'PC-bulk-CTx0m','MC-O-CTx0m',
                    'MC-O-CTx5m','MC-O-CTx8m','MC-O-CTx14m','MC-O-CTx23m'],
                  rotation=45, ha='right')
ax.set_xlabel("")
ax.set_ylabel("Total L1 counts")

plt.savefig("/home/users/wonhee720/Projects/03_CRC/06_results/LINE1/mutrate/TotalL1.PCAWG_NLcolon_added.Chemotherapy_duration.svg")

plt.show()



# All L1s / SNVs in trunk vs. branch in each patient

with open("/home/users/wonhee720/Projects/03_CRC/06_results/SBS_sig_perbranch.pkl", 'rb') as f:
    Branch_SNV = pkl.load(f)

with open("/home/users/wonhee720/Projects/03_CRC/06_results/L1_sig_perbranch.pkl", 'rb') as f:
    Branch_L1 = pkl.load(f)

Branch_SNV_sum = defaultdict(dict)

for pt, dic in Branch_SNV.items():
    Branch_SNV_sum[pt]['trunk'] = []
    Branch_SNV_sum[pt]['branch'] = []
    
    for k, v in dic.items():    
        if pt=='52':
            if "root" in k:
                Branch_SNV_sum[pt]['trunk'].append(v)
            elif "['52-R-FT']" in k:
                Branch_SNV_sum[pt]['trunk'].append(v)
            elif "['52-R2-FT', '52-R-FT']" in k:
                Branch_SNV_sum[pt]['trunk'].append(v)
            elif "['52-R3-FT']" in k: 
                Branch_SNV_sum[pt]['trunk'].append(v)
            elif "['52-R2-FT']" in k:
                Branch_SNV_sum[pt]['trunk'].append(v)
            else:
                Branch_SNV_sum[pt]['branch'].append(v)
                
        else:
            if "root" in k:
                Branch_SNV_sum[pt]['trunk'].append(v)
            else:
                Branch_SNV_sum[pt]['branch'].append(v)


for pt, dic in Branch_SNV_sum.items():
    df = pd.concat(Branch_SNV_sum[pt]['branch']).groupby("Signature").sum("exposure")
    Branch_SNV_sum[pt]['branch'] = df

    df = pd.concat(Branch_SNV_sum[pt]['trunk']).groupby("Signature").sum("exposure")
    Branch_SNV_sum[pt]['trunk'] = df



Branch_L1_sum = defaultdict(dict)

for pt, dic in Branch_L1.items():
    Branch_L1_sum[pt]['trunk'] = []
    Branch_L1_sum[pt]['branch'] = []
    
    for k, v in dic.items():
        if "['31-2019-O6', '31-2019-O5', '31-2019-O7', '31-2019-O4', '31-2019-O3', '31-2019-O2', '31-2019-O1', '31-2019-F', '31-2017-F']" in k:
            Branch_L1_sum[pt]['trunk'].append(v)
        elif "['42-PS-FT', '42-SC-FT', '42-LS6-O7', '42-LS6-O8', '42-LS6-O3', '42-LS6-O6', '42-LS6-O5', '42-LS6-O4', '42-LS6-O2', '42-LS6-O1', '42-LS6-FT']" in k:
            Branch_L1_sum[pt]['trunk'].append(v)
        elif "['43-RT-FT', '43-LLL-O9', '43-LLL-O7', '43-LLL-O4', '43-LLL-O8', '43-LLL-O5', '43-LLL-O3', '43-LLL-FT']" in k:
            Branch_L1_sum[pt]['trunk'].append(v)
        elif "['52-R-FT']" in k:
            Branch_L1_sum[pt]['trunk'].append(v)
        elif "['52-R2-FT', '52-R-FT']" in k:
            Branch_L1_sum[pt]['trunk'].append(v)
        elif "['52-R3-FT']" in k: 
            Branch_L1_sum[pt]['trunk'].append(v)
        elif "['52-R2-FT']" in k:
            Branch_L1_sum[pt]['trunk'].append(v)
        elif "['57-Ovary-O4', '57-Ovary-O3', '57-Ovary-O8', '57-Ovary-O7', '57-Ovary-O2', '57-Ovary-O1', '57-O-FT', '57-LM-O5', '57-LM-O2', '57-LM-O4', '57-LM-O8', '57-LM-O7', '57-LM-O6', '57-LM-O3', '57-LM-O1', '57-CTu-FT', '57-L-FT', '57-CT-O4', '57-CT-O3', '57-CT-O5', '57-CT-O2', '57-CT-O1']" in k:
            Branch_L1_sum[pt]['trunk'].append(v)
        elif "['69-LN-FT', '69-LS78-O5', '69-LS78-O4', '69-LS78-FT', '69-LS4-O3', '69-LS4-O2', '69-LS4-FT', '69-CT-O5', '69-CT-O9', '69-CT-O8', '69-CT-O4', '69-CT-O7', '69-CT-O10', '69-CT-O3', '69-CT-O1', '69-CT-FT']" in k:
            Branch_L1_sum[pt]['trunk'].append(v)
        else:
            Branch_L1_sum[pt]['branch'].append(v)


for pt, dic in Branch_L1_sum.items():
    df = pd.DataFrame(pd.concat(Branch_L1_sum[pt]['branch']).value_counts())
    Branch_L1_sum[pt]['branch'] = df

    df = pd.DataFrame(pd.concat(Branch_L1_sum[pt]['trunk']).value_counts())
    Branch_L1_sum[pt]['trunk'] = df


Total_EPM_L1 = defaultdict(dict)
Total_EPM_L1_df = {}

for pt in ['31','42','43','52','57','69']:
    Total_EPM_L1[pt]['branch'] = {}
    Total_EPM_L1[pt]['trunk'] = {}
    Total_EPM_L1[pt]['branch']['L1'] =  Branch_L1_sum[pt]['branch'].sum()[0]
    Total_EPM_L1[pt]['trunk']['L1'] =  Branch_L1_sum[pt]['trunk'].sum()[0]
    Total_EPM_L1[pt]['branch']['SNV'] =  Branch_SNV_sum[pt]['branch'].loc[Branch_SNV_sum[pt]['branch'].index.isin(['SBS 1',"SBS 5","SBS 40","SBS 3"]) ,:].sum()[0]
    Total_EPM_L1[pt]['trunk']['SNV'] =  Branch_SNV_sum[pt]['trunk'].loc[Branch_SNV_sum[pt]['trunk'].index.isin(['SBS 1',"SBS 5","SBS 40", "SBS 3"]) ,:].sum()[0]


    Total_EPM_L1_df[pt] = ( (Total_EPM_L1[pt]['branch']['L1'] / Total_EPM_L1[pt]['branch']['SNV'])*1000 ) / ( (Total_EPM_L1[pt]['trunk']['L1'] / Total_EPM_L1[pt]['trunk']['SNV'])*1000 ) 


Total_EPM_L1_rate = pd.DataFrame()
for pt, dic in Total_EPM_L1.items():
    trunk = (dic['trunk']['L1']/dic['trunk']['SNV'])*1000
    branch = (dic['branch']['L1']/dic['branch']['SNV'])*1000

    Total_EPM_L1_rate.loc[pt_id_exchanger[pt],'trunk'] = trunk
    Total_EPM_L1_rate.loc[pt_id_exchanger[pt],'branch'] = branch


Total_EPM_L1_rate.sort_index(inplace=True)
Total_EPM_L1_rate.to_csv("/home/users/wonhee720/Projects/03_CRC/06_results/LINE1/mutrate_perbranch/L1_rate.trunk_vs_branch.perpatient.tsv",
                        sep="\t", index=True)


# Average Branch-to-Trunk L1 rate in each patient

L1_per_EPM = pd.DataFrame(Total_EPM_L1_df, index=['AU']).T
fig, ax = plt.subplots(figsize=(9/25.4, 35/25.4) )
sns.barplot(data=L1_per_EPM, y='AU',  ax=ax, capsize=0.3, zorder=5, errorbar='se',
                 err_kws={'color':'black'}, color='#ff7f0e')
#sns.stripplot(data=L1_per_EPM, y='AU', size=2, dodge=True, ax=ax,
#               edgecolor='black', linewidth=0.5, legend=False, zorder=6)

plt.setp(ax.patches, linewidth=0.5, edgecolor='black')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)
ax.spines['left'].set_linewidth(0.5)
ax.spines['bottom'].set_linewidth(0.5)
ax.spines['left'].set_color('black')
ax.spines['bottom'].set_color('black')

plt.ylabel("Branch-to-trunk L1/1000EPMs")

plt.savefig("/home/users/wonhee720/Projects/03_CRC/06_results/LINE1/mutrate/trunk_branch/All.PerPatient.branch-to-trunk.L1perEPMs.svg")
plt.show()








