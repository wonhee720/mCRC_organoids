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


tree_paths = [f"/home/users/wonhee720/Projects/03_CRC/02_vcf/phylogeny/{pt}/outtree_withlength_rescued_allname" for pt in patient_list]


dist_dict_pt = {}
for path, pt in zip(tree_paths, patient_list):
    
    # Read the tree (SNV)
    tree = ete3.Tree(path, format=3)
    
    # Iterate through each nodes except for leaves
    dist_dict = defaultdict(partial(defaultdict, dict))
    dist_dict2 = defaultdict(partial(defaultdict, dict))
    for node in tree.traverse():
        print(f"Node {node.name}!!")
        
        if node.is_leaf():
            pass
        
        else:
            # In each node, find combinations of major x minor descendant
            for i, subnode in enumerate(node.get_children()):
                print(f"Subnode {subnode.name}!!")
                if i ==0:
                    a_desc_count = len(subnode.get_leaves())
                elif i==1:
                    b_desc_count = len(subnode.get_leaves())

                  
                # For each subnode, iterate through its descendant leaves
                for leaf in subnode.get_leaves():
                    dist = subnode.get_distance(target=node, target2=leaf)
                    print(f"Distance from {node.name} to {leaf.name} is {dist}")                   
                    dist_dict[node.name][subnode.name][leaf.name] = dist

                for j, (k, v) in enumerate(dist_dict[node.name].items()):
                    if j==0:
                        first = len(v.items())
                    elif j==1:
                        second = len(v.items())

                # The first is major
                if first >= second :
                    for p, (k, v) in enumerate(dist_dict[node.name].items()):
                        if p==0:
                            dist_dict2[node.name]["Major_lineage"] = dist_dict[node.name][k]
                
                        elif p==1:
                            dist_dict2[node.name]["Minor_lineage"] = dist_dict[node.name][k]
                         
                # The second is major        
                elif first < second:
                    for p, (k, v) in enumerate(dist_dict[node.name].items()):
                        if p==0:
                            dist_dict2[node.name]["Minor_lineage"] = dist_dict[node.name][k]
                            
                        elif p==1:
                            dist_dict2[node.name]["Major_lineage"] = dist_dict[node.name][k]

            # Give descendants ratio info
            if a_desc_count >= b_desc_count:
                dist_dict2[node.name]["Descendants_ratio"] = a_desc_count/b_desc_count
            else:
                dist_dict2[node.name]["Descendants_ratio"] = b_desc_count/a_desc_count
                  
    dist_dict_pt[pt] = dist_dict2


class ComBination:
    def __init__(self, major, minor):
        self.major = major
        self.minor = minor
        
    def scoops(self):
      return list(itertools.product(self.major,self.minor))


asymm_dict_pairwise = defaultdict(dict)
for pt, dist_dict in dist_dict_pt.items():

    for node, v in dist_dict.items():

            major_descendants = v['Major_lineage'].keys()
            minor_descendants = v['Minor_lineage'].keys()

            Combis = ComBination(major_descendants, minor_descendants).scoops()
            # Pairwise combinations (major x minor)
            for (maj,min) in Combis:
                # Access node to leaf distance
                major_dist = dist_dict_pt[pt][node]["Major_lineage"][maj]
                minor_dist = dist_dict_pt[pt][node]["Minor_lineage"][min]

                #1. Mutation difference
                asymm_dict_pairwise[f"{maj}-{min}"]['Asymmetricity diff'] = major_dist-minor_dist

                #2. Mutation ratio 
                try:
                    asymm_dict_pairwise[f"{maj}-{min}"]['Asymmetricity ratio'] = major_dist/minor_dist
                    
                except ZeroDivisionError:
                    print("ZeroDivisionError")
                    asymm_dict_pairwise[f"{maj}-{min}"]['Asymmetricity ratio'] = major_dist

                #3. Descendants number ratio
                asymm_dict_pairwise[f"{maj}-{min}"]['Descendants ratio'] = dist_dict_pt[pt][node]['Descendants_ratio']
            
                if pt == '31':
                    asymm_dict_pairwise[f"{maj}-{min}"]['Chemo duration'] = 14
                    asymm_dict_pairwise[f"{maj}-{min}"]['Pt'] = '31'
                    
                elif pt == '42':
                    asymm_dict_pairwise[f"{maj}-{min}"]['Chemo duration'] = 23
                    asymm_dict_pairwise[f"{maj}-{min}"]['Pt'] = '42'
                    
                elif pt == '57':
                    asymm_dict_pairwise[f"{maj}-{min}"]['Chemo duration'] = 14
                    asymm_dict_pairwise[f"{maj}-{min}"]['Pt'] = '57'
                    
                elif pt == '69':
                    asymm_dict_pairwise[f"{maj}-{min}"]['Chemo duration'] = 5
                    asymm_dict_pairwise[f"{maj}-{min}"]['Pt'] = '69'
                    
                elif pt == '43':
                    asymm_dict_pairwise[f"{maj}-{min}"]['Chemo duration'] = 0
                    asymm_dict_pairwise[f"{maj}-{min}"]['Pt'] = '43'
                    
                elif pt == '52':
                    asymm_dict_pairwise[f"{maj}-{min}"]['Chemo duration'] = 0
                    asymm_dict_pairwise[f"{maj}-{min}"]['Pt'] = '52'
                
                elif pt == '70':
                    asymm_dict_pairwise[f"{maj}-{min}"]['Chemo duration'] = 0
                    asymm_dict_pairwise[f"{maj}-{min}"]['Pt'] = '70'

# Make a dataframe of pairwise comparison for plotting
asymm_pairwise_df = pd.DataFrame(asymm_dict_pairwise,
                                dtype=float).transpose()

asymm_pairwise_df['Descendants ratio'] = asymm_pairwise_df['Descendants ratio'].astype(float)
asymm_pairwise_df['Chemo duration'] = asymm_pairwise_df['Chemo duration'].astype(int)
asymm_pairwise_df['Pt'] = asymm_pairwise_df['Pt'].astype(int)
asymm_pairwise_df['Pt'] = asymm_pairwise_df['Pt'].astype(str)

# Get log value with correponding positivity / negativity
# Fortunately, there is no 0 value in "Asymmetricity diff"
asymm_pairwise_df['Asymmetricity diff_positivity'] = asymm_pairwise_df.apply(lambda row: True if row["Asymmetricity diff"] > 0 else False, axis=1)
asymm_pairwise_df['Asymmetricity diff_log'] = np.log10(np.abs(asymm_pairwise_df['Asymmetricity diff']))
asymm_pairwise_df['Asymmetricity diff_log'] = asymm_pairwise_df.apply(lambda row: row['Asymmetricity diff_log'] if row['Asymmetricity diff_positivity']==True else -row['Asymmetricity diff_log'], axis=1)

# Exclude comparisons including FT, F samples
asymm_pairwise_df = asymm_pairwise_df[~asymm_pairwise_df.index.str.contains("FT|F")]

# Get absolute value of mutation differences for dot size
asymm_pairwise_df['Asymmetricity diff_abs']= np.abs(asymm_pairwise_df['Asymmetricity diff'])

# Patient number assign
asymm_pairwise_df["Pt"] = asymm_pairwise_df["Pt"].replace(pt_exchanger)

# Set manual jitter since scatterplot does not support jitter. If I use stripplot, I can't set sizes of dots according to a specific value.
jitter_strength = 0.15
asymm_pairwise_df['Descendants ratio_jittered'] = asymm_pairwise_df['Descendants ratio'] + np.random.uniform(-jitter_strength, jitter_strength, size=len(asymm_pairwise_df))


fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(100/25.4,120/25.4))

for i, pt in enumerate(["43","52","69","57","31","42"]):
    _df = asymm_pairwise_df_no1_noToT[asymm_pairwise_df_no1_noToT['Pt']==pt_exchanger[pt]]
    
    ymax = np.max(np.abs(_df['Asymmetricity diff_log']))
    
    base_size = (1, 75)
   
    row = i//2
    col = i%2
    scatter = sns.scatterplot(data=_df, x='Descendants ratio',y='Asymmetricity diff_log', legend=False,
               hue="Pt", palette=pt_color_dict_modi,alpha=0.5, s=50,
                 ax=axes[row, col], clip_on=False)

    axes[row,col].set_ylim(-ymax-1,ymax+1)
    if pt=="69":
        axes[row,col].set_xticks([2,2.5,3])
        axes[row,col].set_xticklabels([2,"",3])
    elif pt=="43":
        axes[row,col].set_xticks([2])
        axes[row,col].set_xticklabels([2])    
    
    axes[row,col].axhline(0,linestyle='--', linewidth=0.5, color='black')
    
    if col==0:
        axes[row,col].set_ylabel("Differences in mutation counts\n(Major-Minor) (log)")
    else:
        axes[row,col].set_ylabel("")
    if row==2:
        axes[row,col].set_xlabel("Ratio of descendant clones\n(Major/Minor)")
    else: 
        axes[row,col].set_xlabel("")

    axes[row,col].set_title(pt_exchanger[pt])

plt.tight_layout()
plt.savefig("/home/users/wonhee720/Projects/03_CRC/06_results/asymmetricity/asymmhdiff_log-descenratio(excl_1).excl_difftumor.patient-pairwise.cleanver.svg")
plt.show()  


#### Concordant vs Discordant nodes
fig, ax = plt.subplots(figsize=(40/25.4,40/25.4))

sns.countplot(data=asymm_pairwise_df_no1_noToT, x="Asymmetricity diff_positivity", 
              hue="Asymmetricity diff_positivity", ax=ax, legend=False)

plt.setp(ax.patches, linewidth=0.5, edgecolor='black')
ax.set_xticklabels(labels=["Discordant","Concordant"])
plt.xlabel("")
plt.ylabel("Interclonal comparisons")
ax.axhline(60,linestyle='--', linewidth=0.5, color='black')

plt.savefig("/home/users/wonhee720/Projects/03_CRC/06_results/asymmetricity/asymmdiff-descenratio.excl_difftumor.concordance.pairwise.svg")
plt.show()


chisquare([94, 26], f_exp=[60, 60])



