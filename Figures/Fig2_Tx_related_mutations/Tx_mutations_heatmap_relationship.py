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

AllPt_CTx_df = pd.read_csv("/home/users/wonhee720/Projects/03_CRC/06_results/Tx_related_signature/all_Pt.CTx.tsv", sep="\t", header=0, index_col=0, dtype={"Patient":str})
AllPt_CTx_df_noPt70 = AllPt_CTx_df[~AllPt_CTx_df.index.str.contains("70")]
AllPt_CTx_df_noPt70 = AllPt_CTx_df_noPt70.reindex(id_order_no70)
AllPt_RT_CTx_df_noPt70 = pd.read_csv("/home/users/wonhee720/Projects/03_CRC/06_results/mutcounts.sig_snv_indel_per_sample.from_branch_concat.tsv",
                             sep="\t", index_col=0, header=0)
AllPt_RT_CTx_df_noPt70_subset = pd.read_csv("/home/users/wonhee720/Projects/03_CRC/06_results/mutcounts.sig_snv_indel_per_sample.subset.from_branch_concat.tsv",
                             sep="\t", index_col=0, header=0)


AllPt_RT_CTx_df_noPt70['SBS 5,40'] = AllPt_RT_CTx_df_noPt70['SBS 5'] + AllPt_RT_CTx_df_noPt70['SBS 40']


all_sigs = [ 'ID radiation','SBS cisplatin/oxaliplatin','SBS oxaliplatin','SBS 17b', 'SBS 17a','SBS 25',
  'SBS 5,40', 'SBS 3','SBS 13','SBS 1','SBS 18','SBS 37','SBS 41','Clock-like signature','Chemotherapy-related signature']

AllPt_RT_CTx_df_noPt70['SBS clock-like'] = AllPt_RT_CTx_df_noPt70['SBS 1'] + AllPt_RT_CTx_df_noPt70['SBS 5'] + AllPt_RT_CTx_df_noPt70['SBS 40']

# Added SBS17a to chemotherapy-related signature!!
AllPt_RT_CTx_df_noPt70['Chemotherapy-related signature'] = AllPt_RT_CTx_df_noPt70['Chemotherapy-related mutations']
AllPt_RT_CTx_df_noPt70['Clock-like signature'] = AllPt_RT_CTx_df_noPt70['SBS clock-like']


AllPt_RT_CTx_df_noPt70_subset = AllPt_RT_CTx_df_noPt70[all_sigs].transpose()
sig_color_dict_modi_subet = {}

for x in all_sigs:
    if x == 'SBS oxaliplatin':
        sig_color_dict_modi_subet[x] = sns.color_palette(f"blend:#DBC9C0,{sig_color_dict_modi[x]}", as_cmap=True)
    else:
        sig_color_dict_modi_subet[x] = sns.color_palette(f"blend:#ffffff,{sig_color_dict_modi[x]}", as_cmap=True)
        
vmin = {'SBS 5,40':4000, 'SBS 40':4000, 'SBS 3':490, 'SBS oxaliplatin': 0, # Originally 100 (without fillna 0)
        'SBS cisplatin/oxaliplatin': 0, 'Clock-like signature': 4000, 'Chemotherapy-related signature':200}
vmax = {'SBS 5,40':11500, 'SBS 40':11500, 'SBS 3': 979 ,'SBS oxaliplatin':3000, 'Clock-like signature': 12500, 'Chemotherapy-related signature':4300}
center = {'SBS oxaliplatin':1500} # Originally 400 (without fillna 0)
cbar_yticks = {'ID radiation':850,'SBS cisplatin/oxaliplatin':40,'SBS oxaliplatin':2000, 'SBS 17b':500, 'SBS 17a':500,
  'SBS 5,40':8000,'SBS 40':8000,'SBS 5':900,'SBS 3':750,'SBS 13':700,'SBS 1':2500,'SBS 18':400,'SBS 25':1200,'SBS 37':300,'SBS 41':6000,
               'Clock-like signature':10000,'Chemotherapy-related signature':1500}

AllPt_RT_CTx_df_noPt70_subset.replace(0,np.nan, inplace=True)

# Script making heatmap
sns.set_style("white",  {'figure.facecolor': 'white'})
sigs_num = len(all_sigs)

fig, axes = plt.subplots(sigs_num+2, 2, figsize=(40,15), gridspec_kw={'hspace': 0, 'wspace':0.03, 'width_ratios': [150, 1]})

idx = AllPt_RT_CTx_df_noPt70_subset.index

for i, (index, row) in enumerate(AllPt_RT_CTx_df_noPt70_subset.iterrows()):
    
    if i < (sigs_num-3):
        
        if i==0 :
            sns.heatmap(np.array([AllPt_RT_CTx_df_noPt70_subset.iloc[i,:]]), 
                    ax=axes[i][0], 
                    cmap=sig_color_dict_modi_subet[index], cbar=True, 
                    annot=False, robust=True,
                    cbar_kws = {'drawedges':False, 'shrink':0.1}, cbar_ax = axes[i][1],
                    linewidths = 1, linecolor = 'white',
                    annot_kws={"size": 10})
            
            c_bar = axes[i][0].collections[0].colorbar
            c_bar.set_ticks([cbar_yticks[index]])
            c_bar.set_ticklabels([str(cbar_yticks[index])])

            axes[i][0].set_yticklabels(labels=[idx[i]], rotation=0, fontsize=20)
            
            axes[i][0].set_xticklabels(labels=[id_exchanger[x] for x in AllPt_RT_CTx_df_noPt70_subset.columns], rotation=90, fontsize=20)
            axes[i][0].tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
                
        elif index=='SBS 40' or index=='SBS 3':
            sns.heatmap(np.array([AllPt_RT_CTx_df_noPt70_subset.iloc[i,:]]), 
                xticklabels=False,
                ax=axes[i][0], 
                vmin=vmin[index], vmax=vmax[index],
                cmap=sig_color_dict_modi_subet[index], cbar=True, 
                annot=False, robust=False,
                cbar_kws = {'drawedges':False, 'shrink':0.1}, cbar_ax = axes[i][1],
                linewidths = 1, linecolor = 'white',
                annot_kws={"size": 10})

            c_bar = axes[i][0].collections[0].colorbar
            c_bar.set_ticks([cbar_yticks[index]])
            c_bar.set_ticklabels([str(cbar_yticks[index])])
    
            axes[i][0].set_yticklabels(labels=[idx[i]], rotation=0, fontsize=20)
            
        elif index=='SBS oxaliplatin':
            sns.heatmap(np.array([AllPt_RT_CTx_df_noPt70_subset.iloc[i,:]]), 
                xticklabels=False,
                ax=axes[i][0], 
                vmin=vmin[index], vmax=vmax[index], center=center[index],
                cmap=sig_color_dict_modi_subet[index], cbar=True, 
                annot=False, robust=False,
                cbar_kws = {'drawedges':False, 'shrink':0.1}, cbar_ax = axes[i][1],
                linewidths = 1, linecolor = 'white',
                annot_kws={"size": 10})
            
            c_bar = axes[i][0].collections[0].colorbar
            c_bar.set_ticks([cbar_yticks[index]])
            c_bar.set_ticklabels([str(cbar_yticks[index])])
            
            axes[i][0].set_yticklabels(labels=[idx[i]], rotation=0, fontsize=20)
            
        else:
            sns.heatmap(np.array([AllPt_RT_CTx_df_noPt70_subset.iloc[i,:]]), 
                        xticklabels=False,
                        ax=axes[i][0], 
                        cmap=sig_color_dict_modi_subet[index], cbar=True, 
                        annot=False, robust=True,
                        cbar_kws = {'drawedges':False, 'shrink':0.1}, cbar_ax = axes[i][1],
                        linewidths = 1, linecolor = 'white',
                        annot_kws={"size": 10})

            c_bar = axes[i][0].collections[0].colorbar
            c_bar.set_ticks([cbar_yticks[index]])
            c_bar.set_ticklabels([str(cbar_yticks[index])])
            
            axes[i][0].set_yticklabels(labels=[idx[i]], rotation=0, fontsize=20)
        
    elif (i == (sigs_num-3)):
        sns.heatmap(np.array([AllPt_RT_CTx_df_noPt70_subset.iloc[i,:]]), 
                    xticklabels=False,
                    ax=axes[i][0], 
                    cmap=sig_color_dict_modi_subet[index], cbar=True, 
                    annot=False, robust=True,
                    cbar_kws = {'drawedges':False, 'shrink':0.1}, cbar_ax = axes[i][1],
                    linewidths = 1, linecolor = 'white',
                    annot_kws={"size": 10})
        
        c_bar = axes[i][0].collections[0].colorbar
        c_bar.set_ticks([cbar_yticks[index]])
        c_bar.set_ticklabels([str(cbar_yticks[index])])
        
        
        axes[i][0].set_yticklabels(labels=[idx[i]], rotation=0, fontsize=20)
        
        sns.heatmap(np.array([np.zeros([row.shape[0]])]), ax=axes[i+1][0], cmap=sig_color_dict_modi_subet[index],
                    cbar=False, linewidths = 0,
                    yticklabels=False)        
        
        sns.heatmap(np.array([np.zeros([row.shape[0]])]), ax=axes[i+2][0], cmap=sig_color_dict_modi_subet[index],
                    cbar=False, linewidths = 0,
                    yticklabels=False)        
        
        
    
    else:
        sns.heatmap(np.array([AllPt_RT_CTx_df_noPt70_subset.iloc[i,:]]), 
                    xticklabels=False,
                    ax=axes[i+2][0], 
                    cmap=sig_color_dict_modi_subet[index], cbar=True,
                    vmin=vmin[index], vmax=vmax[index],
                    annot=False, robust=True,
                    cbar_kws = {'drawedges':False, 'shrink':0.1}, cbar_ax = axes[i+2][1],
                    linewidths = 1, linecolor = 'white',
                    annot_kws={"size": 10})
        
        c_bar = axes[i+2][0].collections[0].colorbar
        c_bar.set_ticks([cbar_yticks[index]])
        c_bar.set_ticklabels([str(cbar_yticks[index])])
        
        axes[i+2][0].set_yticklabels(labels=[idx[i]], rotation=0, fontsize=20)
        

for i in np.arange(sigs_num+2):
    if i==0:
            axes[i][0].spines['top'].set_visible(True)
            axes[i][0].spines['right'].set_visible(True)
            axes[i][0].spines['bottom'].set_visible(False)
            axes[i][0].spines['left'].set_visible(True)
            axes[i][0].spines['top'].set_linewidth(3)
            axes[i][0].spines['left'].set_linewidth(3)
            axes[i][0].spines['right'].set_linewidth(3)
            axes[i][0].spines['top'].set_color('black')
            axes[i][0].spines['left'].set_color('black')
            axes[i][0].spines['right'].set_color('black')
            axes[i][0].set_xticklabels(labels=[id_exchanger[x] for x in AllPt_RT_CTx_df_noPt70_subset.columns], rotation=90, fontsize=20)
            axes[i][0].xaxis.set_label_position('top') 
            
    elif (i==(sigs_num-3)):
            axes[i][0].spines['top'].set_visible(False)
            axes[i][0].spines['right'].set_visible(True)
            axes[i][0].spines['bottom'].set_visible(True)
            axes[i][0].spines['left'].set_visible(True)
            axes[i][0].spines['left'].set_linewidth(3)
            axes[i][0].spines['right'].set_linewidth(3)
            axes[i][0].spines['bottom'].set_linewidth(3)
            axes[i][0].spines['left'].set_color('black')
            axes[i][0].spines['right'].set_color('black')
            axes[i][0].spines['bottom'].set_color('black')
            
    elif (i==(sigs_num-2)):
            axes[i][0].spines['top'].set_visible(False)
            axes[i][0].spines['right'].set_visible(False)
            axes[i][0].spines['bottom'].set_visible(False)
            axes[i][0].spines['left'].set_visible(False)

    elif (i==(sigs_num-1)):
            axes[i][0].spines['top'].set_visible(False)
            axes[i][0].spines['right'].set_visible(False)
            axes[i][0].spines['bottom'].set_visible(False)
            axes[i][0].spines['left'].set_visible(False)

            
    elif (i==(sigs_num)):
            axes[i][0].spines['top'].set_visible(True)
            axes[i][0].spines['right'].set_visible(True)
            axes[i][0].spines['bottom'].set_visible(False)
            axes[i][0].spines['left'].set_visible(True)
            axes[i][0].spines['top'].set_linewidth(3)
            axes[i][0].spines['left'].set_linewidth(3)
            axes[i][0].spines['right'].set_linewidth(3)
            axes[i][0].spines['top'].set_color('black')
            axes[i][0].spines['left'].set_color('black')
            axes[i][0].spines['right'].set_color('black')
    
    elif (i==(sigs_num+1)):
            axes[i][0].spines['top'].set_visible(False)
            axes[i][0].spines['right'].set_visible(True)
            axes[i][0].spines['bottom'].set_visible(True)
            axes[i][0].spines['left'].set_visible(True)
            axes[i][0].spines['bottom'].set_linewidth(3)
            axes[i][0].spines['left'].set_linewidth(3)
            axes[i][0].spines['right'].set_linewidth(3)
            axes[i][0].spines['bottom'].set_color('black')
            axes[i][0].spines['left'].set_color('black')
            axes[i][0].spines['right'].set_color('black')
    
    else:
            axes[i][0].spines['top'].set_visible(False)
            axes[i][0].spines['right'].set_visible(True)
            axes[i][0].spines['bottom'].set_visible(False)
            axes[i][0].spines['left'].set_visible(True)
            axes[i][0].spines['left'].set_linewidth(3)
            axes[i][0].spines['right'].set_linewidth(3)
            axes[i][0].spines['left'].set_color('black')
            axes[i][0].spines['right'].set_color('black')
        
plt.savefig("/home/users/wonhee720/Projects/03_CRC/06_results/per_sample/heatmap/SBS_heatmap_allsamples.svg")
plt.show()


# Heatmap for treatment history
Tx_cmap = sns.color_palette(f"blend:#ffffff,#000000", as_cmap=True)

fig, axes = plt.subplots(1, 2, figsize=(40,15/18), gridspec_kw={'hspace': 0, 'wspace':0.03, 'width_ratios': [150, 1]})

sns.heatmap(np.array([AllPt_RT_CTx_df_noPt70['Chemotherapy duration']]),
                ax=axes[0], 
                cmap=Tx_cmap, cbar=True, 
                annot=True, robust=False,
                xticklabels=False,
                cbar_kws = {'drawedges':False, 'shrink':0.1}, cbar_ax = axes[1],
                linewidths = 1, linecolor = 'white',
                annot_kws={"size": 20})

axes[0].set_yticklabels(labels=['Chemotherapy duration (mo)'], rotation=0, fontsize=20)

for spine in ['top', 'right','left','bottom']:
            axes[0].spines[spine].set_visible(True)
            axes[0].spines[spine].set_linewidth(3)
            axes[0].spines[spine].set_color('black')
    
plt.savefig("/home/users/wonhee720/Projects/03_CRC/06_results/per_sample/heatmap/chemoduration_heatmap_allsamples.svg")
plt.show()



# Scatterplot of CTx mutation proportion according to chemotherapy duration per sample

AllPt_CTx_df_noPt_onlyCTx = AllPt_CTx_df_noPt70[AllPt_CTx_df_noPt70['Chemotherapy duration']!=0]

fig, ax = plt.subplots(figsize=(50/25.4, 50/25.4))

#ax = sns.regplot(data=AllPt_CTx_df_noPt70, x='Chemotherapy duration', y='Chemotherapy-related mutations', scatter_kws={'s':0})
ax = sns.swarmplot(data=AllPt_CTx_df_noPt_onlyCTx_noFT, x='Chemotherapy duration', y='Chemotherapy-related mutations', legend=False,
                   hue='Patient', size=5, palette=pt_color_dict)

plt.yticks(ticks=[0,1000,2000,3000,4000,5000,6000], labels=["0","","2,000","","4,000","","6,000"])
    
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)
ax.spines['left'].set_linewidth(0.5)
ax.spines['bottom'].set_linewidth(0.5)
ax.spines['left'].set_color('black')
ax.spines['bottom'].set_color('black')

ax.set_xlabel("Chemotherapy duration (mo)")
ax.set_ylabel("Chemotherapy-related mutations")

#plt.xlim((-0.5,23.5))
#plt.legend(loc='upper left')
plt.tight_layout()

plt.savefig("/home/users/wonhee720/Projects/03_CRC/06_results/Tx_related_signature/Chemo_duration_regression.Chemotherapy_related_mutations.noP435270.noFT.svg")

plt.show()


# Scatterplot of chemo - SBS3 signature
g = sns.lmplot(data=AllPt_CTx_df_Pt31425769_noFT[AllPt_CTx_df_Pt31425769_noFT['Patient']=='69'], x='Chemotherapy-related mutations', height=50/25.4, legend=False,
                   palette = pt_color_dict, y="SBS 3", hue='Patient', scatter_kws={'s':20, 'clip_on':False}, aspect=0.9)

#plt.xlim(0,1500)
plt.xticks(ticks=[0,250,500,750,1000,1250,1500],labels=["0","","500","","1,000","","1,500"],rotation=45, ha='right')
plt.yticks(ticks=[0,250,500,750,1000,1250,1500], labels=["","","500","","1,000","","1,500"])

g.ax.spines['top'].set_visible(False)
g.ax.spines['right'].set_visible(False)
g.ax.spines['bottom'].set_visible(True)
g.ax.spines['left'].set_visible(True)
g.ax.spines['left'].set_linewidth(0.5)
g.ax.spines['bottom'].set_linewidth(0.5)
g.ax.spines['left'].set_color('black')
g.ax.spines['bottom'].set_color('black')
g.ax.set_xlabel("Number of chemotherapy-related mutations")
g.ax.set_ylabel(f"Number of SBS 3 mutations")

plt.savefig(f"/home/users/wonhee720/Projects/03_CRC/06_results/Tx_related_signature/Chemo_mutation.SBS3.CTxPtsOnly.noFT.svg")

plt.show()


pts = ['42','57']

# Scatterplot of chemo - clock-like

for pt in pts:
    
    df = AllPt_CTx_df_Pt31425769_noFT[AllPt_CTx_df_Pt31425769_noFT['Patient']==pt]
    
    g = sns.lmplot(data=df, x='Chemotherapy-related mutations', y='SBS clock-like',
                   hue='Patient', scatter_kws={'s':20, 'clip_on':False }, palette = pt_color_dict, height=50/25.4, legend=False, aspect=0.9 )
    
    if pt == '42':
        #plt.xlim((0,6500))
        plt.xticks(ticks=[0,1500,3000,4500,6000],labels=["0","","3,000","","6,000"],rotation=45, ha='right')
        plt.yticks(ticks=[0,2500,5000,7500,10000,12500], labels=["","2,500","","7,500","","12,500"])
        
    elif pt == '57':
        #plt.xlim((0,1800))
        plt.xticks(ticks=[0,400,800,1200,1600,2000],labels=["0","","800","","1,600",""],rotation=45, ha='right')
        plt.yticks(ticks=[0,2500,5000,7500,10000], labels=["","","5,000","","10,000"])
    
    g.ax.spines['top'].set_visible(False)
    g.ax.spines['right'].set_visible(False)
    g.ax.spines['bottom'].set_visible(True)
    g.ax.spines['left'].set_visible(True)
    g.ax.spines['left'].set_linewidth(0.5)
    g.ax.spines['bottom'].set_linewidth(0.5)
    g.ax.spines['left'].set_color('black')
    g.ax.spines['bottom'].set_color('black')
    
    g.ax.set_xlabel("Number of chemotherapy-related mutations")
    g.ax.set_ylabel("Number of clock-like mutations")
    
    plt.savefig(f"/home/users/wonhee720/Projects/03_CRC/06_results/Tx_related_signature/Chemo_mutation.SBSclock-like.noFT.Pt{pt}.svg")
    
    plt.show()


# Try linear regression
df = AllPt_CTx_df_Pt31425769_noFT[AllPt_CTx_df_Pt31425769_noFT['Patient']=='42']
model = sm.OLS(df['SBS clock-like'],df['Chemotherapy-related mutations'])
results = model.fit()
results.summary()

# Try linear regression
df = AllPt_CTx_df_Pt31425769_noFT[AllPt_CTx_df_Pt31425769_noFT['Patient']=='57']
model = sm.OLS(df['SBS clock-like'],df['Chemotherapy-related mutations'])
results = model.fit()
results.summary()

# Try linear regression
df = AllPt_CTx_df_Pt31425769_noFT[AllPt_CTx_df_Pt31425769_noFT['Patient']=='69']
model = sm.OLS(df['SBS 3'],df['Chemotherapy-related mutations'])
results = model.fit()
results.summary()
