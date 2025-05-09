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
from pprint import pprint

import inhouse_scripts

chromosomes = [str(i) for i in list(range(1, 23)) + ['X', 'Y', 'MT']]
chrom_cat_type = CategoricalDtype(categories=chromosomes, ordered=True)



df = scripts.vcf.tidydf(path)


############# 1st filter ################

rows = df.shape[0]
df['CSQ_CANONICAL'].fillna(".",inplace=True)
df['CSQ_MAX_AF'].fillna(".",inplace=True)
af_ls = []

for i in range(0,rows):
    canonical = df['CSQ_CANONICAL'].iloc[i].split(',')
    index_ls = [n for n,d in enumerate(canonical) if d == "YES"]

    
    if index_ls:
        if "," in str(df['CSQ_MAX_AF'].iloc[i]): 
            af_ls.append(max(df['CSQ_MAX_AF'].iloc[i].split(',')[m] for m in index_ls))
        else:
            af_ls.append(df['CSQ_MAX_AF'].iloc[i])

    else:
        af_ls.append(max(df['CSQ_MAX_AF'].iloc[i].split(",")))
          
df["CSQ_MAX_AF_CANONICAL"] = af_ls
df["CSQ_MAX_AF_CANONICAL"]= df["CSQ_MAX_AF_CANONICAL"].replace(".","0").astype(float)
df['dbSNP'] = df['CSQ_MAX_AF_CANONICAL'] > 0.01

query_string = '(ponSNU_vafPct < 3) and \
(ponBGI_vafPct < 3) and \
(ponPCAWG_vafPct < 3) and \
(ponPCNSL_vafPct < 3)'

df = df.query(query_string, engine='python')
df = df[(df['dbSNP']) & (df[f"vaf:::{normal}"] < 0.03)]



############## 2nd filter ###################


# Filter high overage, low vaf area
q = df[f"total_read_all:::{tumor}"].quantile(q=0.995)
low_vaf = 0.15
low_vafdf = df[df[f"vaf_all:::{tumor}"] < low_vaf]
low_vaf_q = low_vafdf[f"total_read_all:::{tumor}"].quantile(q=0.8)
df = df[~((df[f"total_read_all:::{tumor}"] > q) | 
          ((df[f"vaf_all:::{tumor}"] < low_vaf) & (df[f"total_read_all:::{tumor}"] > low_vaf_q)))]

# Filter with quality (MQ, BQ)
df = df[(df[f"var_read_lowqual:::{tumor}"] + df[f"other_read_all:::{tumor}"] 
       < df[f"var_read_all:::{tumor}"] * 0.75)]

# Filter with MM, Clipfrac
df = df[(df[f"var_meanMM_all:::{tumor}"] < 2.5) & (df[f"var_ClipFrac_all:::{tumor}"] < 0.1)]



########### Singature decomposition ##################

fasta = pysam.FastaFile("/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta")
df[['context_3', 'substitution']] = df.query("VAR_TYPE == 'snp'").apply(lambda x: scripts.sig.context_3(fasta, x['CHROM'], x['POS'], x['ALT']), axis=1, result_type='expand')

result = scripts.sig.sample_decomposition(scripts.sig.data_generation_96(df=df, column_name='context_3', sample_id=k), cancer="Colon", brca=brca_mutation)
sig_result_dfs[pt][k] = result[2]


###### Annotate CCF after 2nd filter ########



############## 3rd filter (selecting clonal mutations only) ###################

tumor = sampleNameBam(f"01_bam_new/{k}.splitmark.realigned.recal.bam")

# Filter out only clonal mutations
df = df[(df['CCF'] != '.') & (df['CCF'] != 'TRUE')]
df = df.astype({"CCF":float})
df  = df[ (df["CCF"] > 0.75) & (df["CCF"] < 1.25)]