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




asymm_dict_pt = {}
for path, pt in zip(tree_paths, patient_list):
    
    # Read the tree (SNV)
    tree = ete3.Tree(path, format=3)
    asymm_dict = defaultdict(dict)
    
    # Iterate through each nodes except for leaves
    for node in tree.traverse():
        print(f"Node {node.name}!!")
        
        if node.is_leaf():
            pass
        
        else:
            dist_avg_ls = []
            descendants_num_ls = []
            # In each node, find children subnodes
            for i, subnode in enumerate(node.get_children()):
                print(f"Subnode {subnode.name}!!")
                
                # Get descendant leaves from the subnode
                descendants_ls = subnode.get_leaves()
                descendants_num = len(descendants_ls)
                descendants_num_ls.append(descendants_num)
                
                dist_sum = 0
                
                # For each subnode, iterate through its descendant leaves
                for leaf in subnode.get_leaves():
                    dist = subnode.get_distance(target=node, target2=leaf)
                    print(f"Distance from {node.name} to {leaf.name} is {dist}")
                    dist_sum += dist
                    
                # Average distance to the leaf from node (in mutations) in the specific subnode 
                dist_avg = dist_sum/(len(list(subnode.get_leaves())))
                print(dist_avg)
                dist_avg_ls.append(dist_avg)
            
            print(dist_avg_ls)
            over = max(dist_avg_ls)
            below = min(dist_avg_ls)
            
            descendants_major = max(descendants_num_ls)
            descendants_minor = min(descendants_num_ls)
            
            asymm_dict[str(node.name)]['Descendants_major'] = descendants_major
            asymm_dict[str(node.name)]['Descendants_minor'] = descendants_minor
            
            # The asymmetricity ratio, difference per node
            asymm_dict[str(node.name)]['Mutations_long'] = over
            asymm_dict[str(node.name)]['Mutations_short'] = below
            asymm_dict[str(node.name)]['Difference'] = over-below
            
            
            try:
                asymm_dict[str(node.name)]['Asymmetricity ratio'] = over/below
            
            except ZeroDivisionError:
                print("ZeroDivisionError")
                asymm_dict[str(node.name)]['Asymmetricity ratio'] = over
            
            if pt == '31':
                asymm_dict[str(node.name)]['Chemo duration'] = 14
                asymm_dict[str(node.name)]['Pt'] = '31'
                
            elif pt == '42':
                asymm_dict[str(node.name)]['Chemo duration'] = 23
                asymm_dict[str(node.name)]['Pt'] = '42'
                
            elif pt == '57':
                asymm_dict[str(node.name)]['Chemo duration'] = 14
                asymm_dict[str(node.name)]['Pt'] = '57'
                
            elif pt == '69':
                asymm_dict[str(node.name)]['Chemo duration'] = 5
                asymm_dict[str(node.name)]['Pt'] = '69'
                
            elif pt == '43':
                asymm_dict[str(node.name)]['Chemo duration'] = 0
                asymm_dict[str(node.name)]['Pt'] = '43'
                
            elif pt == '52':
                asymm_dict[str(node.name)]['Chemo duration'] = 0
                asymm_dict[str(node.name)]['Pt'] = '52'
            
            elif pt == '70':
                asymm_dict[str(node.name)]['Chemo duration'] = 0
                asymm_dict[str(node.name)]['Pt'] = '70'
            
    asymm_dict_pt[pt] = asymm_dict


# Tier system => suitable for per patient analysis
asymm_dict_pt['31']['31-2019-O6_ 31-2019-O5_ 31-2019-O7_ 31-2019-O4_ 31-2019-O3_ 31-2019-O2_ 31-2019-O1_ 31-2019-F_ 31-2017-F']['tier'] = 1
asymm_dict_pt['31']['31-2019-O6_ 31-2019-O5_ 31-2019-O7_ 31-2019-O4_ 31-2019-O3_ 31-2019-O2_ 31-2019-O1_ 31-2019-F']['tier'] = 2
asymm_dict_pt['31']['31-2019-O6_ 31-2019-O5_ 31-2019-O7_ 31-2019-O4_ 31-2019-O3_ 31-2019-O2_ 31-2019-O1']['tier'] = 3
asymm_dict_pt['31']['31-2019-O5_ 31-2019-O7_ 31-2019-O4_ 31-2019-O3_ 31-2019-O2_ 31-2019-O1']['tier'] = 4
asymm_dict_pt['31']['31-2019-O5_ 31-2019-O7_ 31-2019-O4_ 31-2019-O3_ 31-2019-O2']['tier'] = 5
asymm_dict_pt['31']['31-2019-O7_ 31-2019-O4_ 31-2019-O3_ 31-2019-O2']['tier'] = 6
asymm_dict_pt['31']['31-2019-O7_ 31-2019-O4_ 31-2019-O3']['tier'] = 7
asymm_dict_pt['31']['31-2019-O4_ 31-2019-O3']['tier'] = 8

asymm_dict_pt['42']['42-PS-FT_ 42-SC-FT_ 42-LS6-O7_ 42-LS6-O8_ 42-LS6-O3_ 42-LS6-O6_ 42-LS6-O5_ 42-LS6-O4_ 42-LS6-O2_ 42-LS6-O1_ 42-LS6-FT']['tier'] = 1
asymm_dict_pt['42']['42-SC-FT_ 42-LS6-O7_ 42-LS6-O8_ 42-LS6-O3_ 42-LS6-O6_ 42-LS6-O5_ 42-LS6-O4_ 42-LS6-O2_ 42-LS6-O1_ 42-LS6-FT']['tier'] = 2
asymm_dict_pt['42']['42-LS6-O7_ 42-LS6-O8_ 42-LS6-O3_ 42-LS6-O6_ 42-LS6-O5_ 42-LS6-O4_ 42-LS6-O2_ 42-LS6-O1_ 42-LS6-FT']['tier'] = 3
asymm_dict_pt['42']['42-LS6-O7_ 42-LS6-O8_ 42-LS6-O3_ 42-LS6-O6_ 42-LS6-O5_ 42-LS6-O4_ 42-LS6-O2_ 42-LS6-O1']['tier'] = 4
asymm_dict_pt['42']['42-LS6-O7_ 42-LS6-O8_ 42-LS6-O3_ 42-LS6-O6_ 42-LS6-O5_ 42-LS6-O4_ 42-LS6-O2']['tier'] = 5
asymm_dict_pt['42']['42-LS6-O8_ 42-LS6-O3_ 42-LS6-O6_ 42-LS6-O5_ 42-LS6-O4_ 42-LS6-O2']['tier'] = 6
asymm_dict_pt['42']['42-LS6-O3_ 42-LS6-O6_ 42-LS6-O5_ 42-LS6-O4_ 42-LS6-O2']['tier'] = 7
asymm_dict_pt['42']['42-LS6-O6_ 42-LS6-O5_ 42-LS6-O4_ 42-LS6-O2']['tier'] = 8
asymm_dict_pt['42']['42-LS6-O5_ 42-LS6-O4_ 42-LS6-O2']['tier'] = 9
asymm_dict_pt['42']['42-LS6-O4_ 42-LS6-O2']['tier'] = 10

asymm_dict_pt['43']['43-RT-FT_ 43-LLL-O9_ 43-LLL-O7_ 43-LLL-O4_ 43-LLL-O8_ 43-LLL-O5_ 43-LLL-O3_ 43-LLL-FT']['tier'] = 1
asymm_dict_pt['43']['43-LLL-O9_ 43-LLL-O7_ 43-LLL-O4_ 43-LLL-O8_ 43-LLL-O5_ 43-LLL-O3_ 43-LLL-FT']['tier'] = 2
asymm_dict_pt['43']['43-LLL-O9_ 43-LLL-O7_ 43-LLL-O4_ 43-LLL-O8_ 43-LLL-O5_ 43-LLL-O3']['tier'] = 3
asymm_dict_pt['43']['43-LLL-O9_ 43-LLL-O7_ 43-LLL-O4']['tier'] = 4
asymm_dict_pt['43']['43-LLL-O8_ 43-LLL-O5_ 43-LLL-O3']['tier'] = 4
asymm_dict_pt['43']['43-LLL-O7_ 43-LLL-O4']['tier'] = 5
asymm_dict_pt['43']['43-LLL-O8_ 43-LLL-O5']['tier'] = 5

asymm_dict_pt['52']['52-R2-FT_ 52-R-FT_ 52-R3-FT_ 52-LRt-O7_ 52-LRt-O5_ 52-LRt-O2_ 52-LRt-O6_ 52-LRt-O4_ 52-LRt-O1_ 52-LRT-FT']['tier'] = 1
asymm_dict_pt['52']['52-R2-FT_ 52-R-FT']['tier'] = 2
asymm_dict_pt['52']['52-R3-FT_ 52-LRt-O7_ 52-LRt-O5_ 52-LRt-O2_ 52-LRt-O6_ 52-LRt-O4_ 52-LRt-O1_ 52-LRT-FT']['tier'] = 2
asymm_dict_pt['52']['52-LRt-O7_ 52-LRt-O5_ 52-LRt-O2_ 52-LRt-O6_ 52-LRt-O4_ 52-LRt-O1_ 52-LRT-FT']['tier'] = 3
asymm_dict_pt['52']['52-LRt-O7_ 52-LRt-O5_ 52-LRt-O2_ 52-LRt-O6_ 52-LRt-O4_ 52-LRt-O1']['tier'] = 4
asymm_dict_pt['52']['52-LRt-O5_ 52-LRt-O2_ 52-LRt-O6_ 52-LRt-O4_ 52-LRt-O1']['tier'] = 5
asymm_dict_pt['52']['52-LRt-O5_ 52-LRt-O2']['tier'] = 6
asymm_dict_pt['52']['52-LRt-O6_ 52-LRt-O4_ 52-LRt-O1']['tier'] = 6
asymm_dict_pt['52']['52-LRt-O4_ 52-LRt-O1']['tier'] = 7

asymm_dict_pt['57']['57-Ovary-O4_ 57-Ovary-O3_ 57-Ovary-O8_ 57-Ovary-O7_ 57-Ovary-O2_ 57-Ovary-O1_ 57-O-FT_ 57-LM-O5_ 57-LM-O2_ 57-LM-O4_ 57-LM-O8_ 57-LM-O7_ 57-LM-O6_ 57-LM-O3_ 57-LM-O1_ 57-CTu-FT_ 57-L-FT_ 57-CT-O4_ 57-CT-O3_ 57-CT-O5_ 57-CT-O2_ 57-CT-O1']['tier'] = 1
asymm_dict_pt['57']['57-Ovary-O4_ 57-Ovary-O3_ 57-Ovary-O8_ 57-Ovary-O7_ 57-Ovary-O2_ 57-Ovary-O1_ 57-O-FT']['tier'] = 2
asymm_dict_pt['57']['57-LM-O5_ 57-LM-O2_ 57-LM-O4_ 57-LM-O8_ 57-LM-O7_ 57-LM-O6_ 57-LM-O3_ 57-LM-O1_ 57-CTu-FT_ 57-L-FT_ 57-CT-O4_ 57-CT-O3_ 57-CT-O5_ 57-CT-O2_ 57-CT-O1']['tier'] = 2
asymm_dict_pt['57']['57-Ovary-O4_ 57-Ovary-O3_ 57-Ovary-O8_ 57-Ovary-O7_ 57-Ovary-O2_ 57-Ovary-O1']['tier'] = 3
asymm_dict_pt['57']['57-LM-O5_ 57-LM-O2_ 57-LM-O4_ 57-LM-O8_ 57-LM-O7_ 57-LM-O6_ 57-LM-O3_ 57-LM-O1_ 57-CTu-FT']['tier'] = 3
asymm_dict_pt['57']['57-L-FT_ 57-CT-O4_ 57-CT-O3_ 57-CT-O5_ 57-CT-O2_ 57-CT-O1']['tier'] = 3
asymm_dict_pt['57']['57-Ovary-O3_ 57-Ovary-O8_ 57-Ovary-O7_ 57-Ovary-O2_ 57-Ovary-O1']['tier'] = 4
asymm_dict_pt['57']['57-LM-O5_ 57-LM-O2_ 57-LM-O4_ 57-LM-O8_ 57-LM-O7_ 57-LM-O6_ 57-LM-O3_ 57-LM-O1']['tier'] = 4
asymm_dict_pt['57']['57-CT-O4_ 57-CT-O3_ 57-CT-O5_ 57-CT-O2_ 57-CT-O1']['tier'] = 4
asymm_dict_pt['57']['57-Ovary-O3_ 57-Ovary-O8_ 57-Ovary-O7_ 57-Ovary-O2']['tier'] = 5
asymm_dict_pt['57']['57-LM-O2_ 57-LM-O4_ 57-LM-O8_ 57-LM-O7_ 57-LM-O6_ 57-LM-O3_ 57-LM-O1']['tier'] = 5
asymm_dict_pt['57']['57-CT-O3_ 57-CT-O5_ 57-CT-O2_ 57-CT-O1']['tier'] = 5
asymm_dict_pt['57']['57-Ovary-O8_ 57-Ovary-O7_ 57-Ovary-O2']['tier'] = 6
asymm_dict_pt['57']['57-LM-O4_ 57-LM-O8_ 57-LM-O7_ 57-LM-O6_ 57-LM-O3_ 57-LM-O1']['tier'] = 6
asymm_dict_pt['57']['57-CT-O3_ 57-CT-O5_ 57-CT-O2']['tier'] = 6
asymm_dict_pt['57']['57-Ovary-O7_ 57-Ovary-O2']['tier'] = 7
asymm_dict_pt['57']['57-LM-O8_ 57-LM-O7_ 57-LM-O6_ 57-LM-O3_ 57-LM-O1']['tier'] = 7
asymm_dict_pt['57']['57-CT-O5_ 57-CT-O2']['tier'] = 7
asymm_dict_pt['57']['57-LM-O8_ 57-LM-O7_ 57-LM-O6_ 57-LM-O3']['tier'] = 8
asymm_dict_pt['57']['57-LM-O7_ 57-LM-O6_ 57-LM-O3']['tier'] = 9
asymm_dict_pt['57']['57-LM-O7_ 57-LM-O6']['tier'] = 10

asymm_dict_pt['69']['69-LN-FT_ 69-LS78-O5_ 69-LS78-O4_ 69-LS78-FT_ 69-LS4-O3_ 69-LS4-O2_ 69-LS4-FT_ 69-CT-O5_ 69-CT-O9_ 69-CT-O8_ 69-CT-O4_ 69-CT-O7_ 69-CT-O10_ 69-CT-O3_ 69-CT-O1_ 69-CT-FT']['tier'] = 1
asymm_dict_pt['69']['69-LN-FT_ 69-LS78-O5_ 69-LS78-O4_ 69-LS78-FT_ 69-LS4-O3_ 69-LS4-O2_ 69-LS4-FT_ 69-CT-O5_ 69-CT-O9_ 69-CT-O8_ 69-CT-O4_ 69-CT-O7_ 69-CT-O10_ 69-CT-O3_ 69-CT-O1']['tier'] = 2
asymm_dict_pt['69']['69-LS78-O5_ 69-LS78-O4_ 69-LS78-FT_ 69-LS4-O3_ 69-LS4-O2_ 69-LS4-FT_ 69-CT-O5_ 69-CT-O9_ 69-CT-O8_ 69-CT-O4_ 69-CT-O7_ 69-CT-O10_ 69-CT-O3_ 69-CT-O1']['tier'] = 3
asymm_dict_pt['69']['69-LS78-O5_ 69-LS78-O4_ 69-LS78-FT_ 69-LS4-O3_ 69-LS4-O2_ 69-LS4-FT']['tier'] = 4
asymm_dict_pt['69']['69-CT-O5_ 69-CT-O9_ 69-CT-O8_ 69-CT-O4_ 69-CT-O7_ 69-CT-O10_ 69-CT-O3_ 69-CT-O1']['tier'] = 4
asymm_dict_pt['69']['69-LS78-O5_ 69-LS78-O4_ 69-LS78-FT']['tier'] = 5
asymm_dict_pt['69']['69-LS4-O3_ 69-LS4-O2_ 69-LS4-FT']['tier'] = 5
asymm_dict_pt['69']['69-CT-O5_ 69-CT-O9_ 69-CT-O8_ 69-CT-O4']['tier'] = 5
asymm_dict_pt['69']['69-CT-O7_ 69-CT-O10_ 69-CT-O3_ 69-CT-O1']['tier'] = 5
asymm_dict_pt['69']['69-LS78-O5_ 69-LS78-O4']['tier'] = 6
asymm_dict_pt['69']['69-LS4-O2_ 69-LS4-FT']['tier'] = 6
asymm_dict_pt['69']['69-CT-O9_ 69-CT-O8_ 69-CT-O4']['tier'] = 6
asymm_dict_pt['69']['69-CT-O7_ 69-CT-O10']['tier'] = 6
asymm_dict_pt['69']['69-CT-O3_ 69-CT-O1']['tier'] = 6
asymm_dict_pt['69']['69-CT-O8_ 69-CT-O4']['tier'] = 7

asymm_dict_pt['70']['70-ColonT-FT_ 70-Colon-T_O6_ 70-Colon-T_O3_ 70-Colon-T_O9_ 70-Colon-T_O5_ 70-Colon-T_O10_ 70-Colon-T_O7_ 70-Colon-T_O4_ 70-Colon-T_O2_ 70-Colon-T_O8_ 70-Colon-T_O1']['tier'] = 1
asymm_dict_pt['70']['70-Colon-T_O6_ 70-Colon-T_O3_ 70-Colon-T_O9_ 70-Colon-T_O5_ 70-Colon-T_O10_ 70-Colon-T_O7_ 70-Colon-T_O4_ 70-Colon-T_O2_ 70-Colon-T_O8_ 70-Colon-T_O1']['tier'] = 2
asymm_dict_pt['70']['70-Colon-T_O3_ 70-Colon-T_O9_ 70-Colon-T_O5_ 70-Colon-T_O10_ 70-Colon-T_O7_ 70-Colon-T_O4_ 70-Colon-T_O2_ 70-Colon-T_O8_ 70-Colon-T_O1']['tier'] = 3
asymm_dict_pt['70']['70-Colon-T_O9_ 70-Colon-T_O5_ 70-Colon-T_O10_ 70-Colon-T_O7_ 70-Colon-T_O4_ 70-Colon-T_O2_ 70-Colon-T_O8_ 70-Colon-T_O1']['tier'] = 4
asymm_dict_pt['70']['70-Colon-T_O5_ 70-Colon-T_O10_ 70-Colon-T_O7_ 70-Colon-T_O4_ 70-Colon-T_O2_ 70-Colon-T_O8_ 70-Colon-T_O1']['tier'] = 5
asymm_dict_pt['70']['70-Colon-T_O5_ 70-Colon-T_O10']['tier'] = 6
asymm_dict_pt['70']['70-Colon-T_O7_ 70-Colon-T_O4_ 70-Colon-T_O2_ 70-Colon-T_O8_ 70-Colon-T_O1']['tier'] = 6
asymm_dict_pt['70']['70-Colon-T_O7_ 70-Colon-T_O4_ 70-Colon-T_O2']['tier'] = 7
asymm_dict_pt['70']['70-Colon-T_O8_ 70-Colon-T_O1']['tier'] = 7
asymm_dict_pt['70']['70-Colon-T_O7_ 70-Colon-T_O4']['tier'] = 8



for i, v in enumerate(asymm_dict_pt.values()):
    if i==0:
        asymm_dict_node=copy.deepcopy(v)
    else:
        asymm_dict_node.update(v)



# Make a dataframe for plotting
asymm_df = pd.DataFrame(asymm_dict_node)
asymm_df = asymm_df.transpose()
asymm_df[['Chemo duration','Asymmetricity ratio','tier','Difference']] = asymm_df[['Chemo duration','Asymmetricity ratio','tier','Difference']].astype(float) 
asymm_df['Asymmetricity ratio_log'] = np.log10(asymm_df['Asymmetricity ratio'])
asymm_df['Difference_log'] = np.log10(asymm_df['Difference'])
asymm_df['Chemotherapy'] = asymm_df['Chemo duration']>0


asymm_df_CTx = asymm_df[asymm_df['Chemotherapy']==True]
asymm_df_noCTx = asymm_df[asymm_df['Chemotherapy']==False]

asymm_df_pt = {}

for pt in patient_id_list:
    asymm_df_pt[pt] = asymm_df[asymm_df['Pt']==pt]



# Dictionary of average branch mutation counts per patient 
each_branch_mutnum_dict = {}
for pt in patient_list:
    if pt == '52':
        # Get sum except for root in each primaries for pt 52
        for i, (k, v) in enumerate(rescued_brs_dfs[pt].items()):        
            if k == "52.root":
                pass
            elif k == "['52-R2-FT']":
                pass
            elif k == "['52-R-FT']":
                pass
            elif k == "['52-R3-FT', '52-LRt-O7', '52-LRt-O5', '52-LRt-O2', '52-LRt-O6', '52-LRt-O4', '52-LRt-O1', '52-LRT-FT']":
                pass
            else:
                branch_mut_num = v.shape[0]
                each_branch_mutnum_dict[k] = branch_mut_num
        
    else:    
        # Get sum except for root
        for i, (k, v) in enumerate(rescued_brs_dfs[pt].items()):        
            ptroot = pt+'.root'
            if k == ptroot:
                pass
            else:
                branch_mut_num = v.shape[0]
                each_branch_mutnum_dict[k] = branch_mut_num



# Pie graph for asymmetricity under chemotherapy in organoid branches

pie_df = asymm_df[asymm_df['Pt'].str.contains('42|57|69')][['Asymmetricity ratio','Mutations_long','Mutations_short','Descendants_major','Descendants_minor']]

for index, row in pie_df.iterrows():

    descendants = [int(row['Descendants_minor']), int(row['Descendants_major'])]
    ###### Have to check because not all major lineages have long mutations #######
    ###### Have to check places with 0 or 999 values #######
    tmp_dict = {'minor':{'desc':int(row['Descendants_minor']), 'muts':int(row['Mutations_short'])}, 
                     'major':{'desc':int(row['Descendants_major']), 'muts':int(row['Mutations_long'])} }

    major_radius = float(row['Asymmetricity ratio'])
    
    # Limit transparency for too small branches
    if major_radius > 5:
        major_radius = 5
    
    fig, ax = plt.subplots()
    
    wedges_major,label_major,text_major  = ax.pie(descendants, radius=1, labels=['Minor', 'Major'], 
                                                  startangle=90, counterclock=False, colors=['#ED1C24','#009245'], 
                                                  wedgeprops = {'linewidth': 3, 'edgecolor':'black'},
                                                  autopct=str(tmp_dict['major']['desc']) +"\n"+ "(" + str(tmp_dict['major']['muts']) + ")",
                                                  pctdistance = 0.7, 
                                                  textprops = {'size':30}, labeldistance=None)
    wedges_major[0].set_visible(False)
    #label_major[0].set_visible(False)
    text_major[0].set_visible(False)
    
    wedges_minor,label_minor,text_minor = ax.pie(descendants, radius=1, labels=['Minor', 'Major'], explode=[0.4,0] ,
                                                 startangle=90, counterclock=False, colors=['#ED1C24','#009245'], 
                                                 wedgeprops = {'linewidth': 3, 'edgecolor':'black','alpha':1/major_radius}, 
                                                 autopct=str(tmp_dict['minor']['desc']) +"\n"+ "(" + str(tmp_dict['minor']['muts']) + ")",
                                                 pctdistance = 0.7, 
                                                 textprops = {'size':30}, labeldistance=None)
    wedges_minor[1].set_visible(False)
    #label_minor[1].set_visible(False)
    text_minor[1].set_visible(False)
    
    plt.savefig(f"/home/users/wonhee720/Projects/03_CRC/06_results/asymmetricity/asymm_piegraph/{row.name}.svg")
    print(row.name)
    plt.show()




















