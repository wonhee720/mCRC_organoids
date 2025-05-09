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
from scipy.stats import ttest_ind
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
import ete3
import math
import matplotlib.font_manager as fm
from matplotlib.ticker import MaxNLocator
import matplotlib
import statsmodels.api as sm
import statistics
from scipy import stats
from scipy.stats import f_oneway
from scipy.stats import chisquare
import scipy.stats as stat
from matplotlib.tri import Triangulation
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import itertools
import colorcet as cc
import logging

import inhouse_scripts

chromosomes = [str(i) for i in list(range(1, 23)) + ['X', 'Y']]
chrom_cat_type = CategoricalDtype(categories=chromosomes, ordered=True)


LINE1_df_final = pd.read_csv("/home/users/wonhee720/Projects/03_CRC/02_vcf/LINE1/review_final_all_samples_repos_correct_whole.tumor.v12.txt",
                  sep="\t",header=0) # Dataframe from manual review of L1 calls


LINE1_df_final_soL1 = LINE1_df_final[LINE1_df_final['L1_type']=='soL1R']


LINE1_df_final_soL1['L1_info_5'] = LINE1_df_final_soL1['L1_info'].str.split("_", expand=True)[1]
LINE1_df_final_soL1['L1_info_3'] = LINE1_df_final_soL1['L1_info'].str.split("_", expand=True)[2]

soL1R_bpts = pd.DataFrame(LINE1_df_final_soL1.value_counts(['CHROM','bkp1','bkp2'])).reset_index()

tmp_df_ls = []

for i in range(soL1R_bpts.shape[0]):
    chrom = soL1R_bpts.loc[i,'CHROM']
    bpt1 = soL1R_bpts.loc[i,'bkp1']
    bpt2 = soL1R_bpts.loc[i,'bkp2']

    tmp_df = LINE1_df_final_soL1[(LINE1_df_final_soL1['CHROM']==chrom) & (LINE1_df_final_soL1['bkp1']==bpt1) & (LINE1_df_final_soL1['bkp2']==bpt2)]

    L1size_min = float(tmp_df['L1_size'].min())
    L1size_max = float(tmp_df['L1_size'].max())

    if L1size_max-L1size_min > 1500 :
        tmp_df['L1_size_discrepant'] = True

    else:
        tmp_df[['L1_info_5','L1_info_3']] = tmp_df[['L1_info_5','L1_info_3']].astype(float)
        L15 = np.nanmin(tmp_df['L1_info_5'])
        L13 = np.nanmax(tmp_df['L1_info_3'])
        tmp_df['L1_info'] = "L1C" + "_" + str(L15) + "_" + str(L13)
        tmp_df['L1_size'] = L13-L15
        tmp_df['L1_size_discrepant'] = False

    tmp_df_ls.append(tmp_df)

soL1R_fixed = pd.concat(tmp_df_ls)

soL1R_fixed = soL1R_fixed.sort_values(['CHROM','bkp1','bkp2'])

LINE1_df_final_rest = LINE1_df_final[LINE1_df_final['L1_type']!='soL1R']

LINE1_df_final = pd.concat([soL1R_fixed,LINE1_df_final_rest]).sort_values(['CHROM','bkp1','bkp2'])


LINE1_df_final.to_csv("/home/users/wonhee720/Projects/03_CRC/02_vcf/LINE1/review_final_all_samples_repos_correct_whole.tumor.v13.tsv", 
                   sep="\t", index=False) # Save the final dataframe

###### AFter the final dataframe save, additional manual review done to rescue unfit L1 calls #######












































