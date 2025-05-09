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

snv_filtered_3rd_tsv_path_dict = {k : [f"/home/users/wonhee720/Projects/03_CRC/02_vcf/union/{x}/{x}.snps.filtered_3rd.seqzcn-ccf.tsv" for x in sample_id_dict[k]] for k in patient_list }
snv_filtered_3rd_vcf_path_dict = {k : [f"/home/users/wonhee720/Projects/03_CRC/02_vcf/union/{x}/{x}.snps.filtered_3rd.seqzcn-ccf.vcf.gz" for x in sample_id_dict[k]] for k in patient_list }

dfs_snv_filtered_3rd = defaultdict(dict)

for k,v in snv_filtered_3rd_vcf_path_dict.items():
    print("patient " + str(k) + " !!!")
    dfs_snv_filtered_3rd[k] = {}
    for sample, path in zip(sample_id_dict[k], v):
        dfs_snv_filtered_3rd[k][sample] = scripts.vcf.tidydf(path)
        print(path + " loaded !!!")


# Making merged dataframe

for patient in dfs_snv_filtered_3rd.keys():
    
    # Lines for merging clonal organoids only
    #sample_id_list_clonal=[]
    #for sample in sample_id_dict[patient]:
    #    if "_O" in str(sample) or "42" in str(sample):
    #        sample_id_list_clonal.append(sample)
            
    #dfs_snv_filtered_2nd_clonal = {key: dfs_snv_filtered_2nd[patient][key] for key in sample_id_list_clonal}
    
    sample_id_list = sample_id_dict[patient]
    
    for i, (k, v) in enumerate(dfs_snv_filtered_3rd[patient].items()):
        df = v.df
        tumor = sampleNameBam(f"01_bam_new/{k}.splitmark.realigned.recal.bam")
        df[k] = df[f"vaf_all:::{tumor}"]
        if i==0:
            merge_df = df[["CHROM","POS","REF","ALT",k]]
            
        else:
            merge_df = merge_df.merge(df[["CHROM","POS","REF","ALT", k]], how="outer", on=["CHROM","POS","REF","ALT"])
    
    ls = sample_id_list+["CHROM","POS","REF","ALT"]
    merge_df = merge_df[ls]

    merge_df["ID"] = '.'
    merge_df["QUAL"] = 0
    merge_df["FILTER"] = "PASS"
    merge_df["AS_FilterStatus"] = "SITE"
    
    merge_df = merge_df[["CHROM","POS","ID","REF","ALT","QUAL","FILTER","AS_FilterStatus"]]
    merge_df.rename(columns={"CHROM":"#CHROM"}, inplace=True)
    merge_df.sort_values(by=["#CHROM","POS"], inplace=True)
    print(patient)
    display(merge_df.head())
    print(merge_df.shape)
    merge_df.to_csv(f"/home/users/wonhee720/Projects/03_CRC/02_vcf/phylogeny/{patient}/{patient}.snps.clonal.merged.tsv", sep="\t",index=False)



####### Annotate the merged dataframe with bam information from each sample ######


snv_merged_ori ={}
for pt in patient_list:
    snv_merged_ori[pt] = scripts.vcf.tidydf(f"/home/users/wonhee720/Projects/03_CRC/02_vcf/phylogeny/{pt}/{pt}.snps.clonal.merged.readinfo.vcf.gz").df


vaf_cutoff_dict = {'31':0.025, '42': 0.025, '43':0.025,'52': 0.025, '57': 0.025, '69':0.025}


# making sorted index (by sample order) for heatmap visualization
merged_for_index = {}
for pt in patient_list:
    _dfs = dfs_snv_filtered_3rd[pt]
    
    for i, (k, v) in enumerate(_dfs.items()):
        df = v.df
        tumor = sampleNameBam(f"01_bam_new/{k}.splitmark.realigned.recal.bam")
        df[k] = df[f"vaf_all:::{tumor}"]
        if i==0:
            merged_for_index[pt] = df[["CHROM","POS","REF","ALT",k]]
            
        else:
            merged_for_index[pt] = merged_for_index[pt].merge(df[["CHROM","POS","REF","ALT",k]], how="outer", on=["CHROM","POS","REF","ALT"])
    
    sample_id_list = list(_dfs.keys())
    ls = ["CHROM","POS","REF","ALT"] + sample_id_list
    merged_for_index[pt] = merged_for_index[pt][ls]
    
    merged_for_index[pt].set_index(keys=["CHROM","POS","REF","ALT"], inplace=True)


# Dictionary for mutations excluding commonly shared ones per patient
snv_merged = copy.deepcopy(snv_merged_ori)

# Dictionary for commonly sharing mutations per patient
snv_root = {}

for pt in patient_list:
    _df = snv_merged[pt]
    total_count = _df.shape
    print(f"The total number of merged SNVs in patient {pt} is {total_count} !!!")
    
    # Line for including clonal organoid samples only
    #sample_id_list_clonal=[]
    #for sample in sample_id_dict[pt]:
    #    if "_O" in str(sample) or "42" in str(sample):
    #        sample_id_list_clonal.append(sample)
            
    # Filter mutations with total (var+ref) read count of 5
    vaf_name_list = []
    for k in sample_id_dict[pt]:
        tumor = sampleNameBam(f"01_bam_new/{k}.splitmark.realigned.recal.bam")
        _df[f"sum_read_all:::{tumor}"] = _df[f"var_read_all:::{tumor}"] + _df[f"ref_read_all:::{tumor}"]
        _df = _df[_df[f"sum_read_all:::{tumor}"] > 5]
        vaf_name_list.append(f'vaf_all:::{tumor}')
    
    ls = ["CHROM","POS","REF","ALT"] + vaf_name_list
    _df = _df[ls]
    
    _df.set_index(keys=["CHROM","POS","REF","ALT"], inplace=True)
    
    _df = _df[vaf_name_list]
    
    
    # Make matrix of 0 and 1 using VAF cutoff designated for each patient
    for k, vaf_name in zip(sample_id_dict[pt], vaf_name_list):
        
        # Set special cutoffs
        #if k == "42-PS-FT": 
        #    _df[k] = (_df[vaf_name] > 0.02)
        #    
        #else:
        #    _df[k] = (_df[vaf_name] > vaf_cutoff_dict[pt])
        
        _df[k] = (_df[vaf_name] > vaf_cutoff_dict[pt])
    
    counts=len(sample_id_dict[pt])

    _df["sum_vaf0"] = _df[sample_id_dict[pt]].sum(axis=1)

    
    # Heatmap before excluding mutations shared in all samples
    fig, ax = plt.subplots(figsize=(8,15))
    fig.suptitle(f"Heatmap of number {pt} patient")
    sns.heatmap(data = _df[sample_id_dict[pt]], cbar=False, cmap=["white","steelblue"], ax=ax)
    plt.show()    
    
    
    # Exclude mutations shared in all samples
    filtered_snv_merged_df = _df[_df["sum_vaf0"]<counts] 
    snv_root[pt] = _df[_df["sum_vaf0"]==counts]
    
    # VAF plot excluding mutations shared in all samples
    for vaf_name in vaf_name_list:
        fig, ax = plt.subplots()
        sns.distplot(filtered_snv_merged_df[vaf_name], bins=100)
        ax.set_xlim((0,1.2))
        ax.set_ylim((0,10))
        plt.show()
    
    # Show number of mutations in the matrix
    display(filtered_snv_merged_df["sum_vaf0"].value_counts())
    
    filtered_snv_merged_df = filtered_snv_merged_df.iloc[:,counts:-1]
    
    filtered_snv_merged_df = merged_for_index[pt][[]].merge(filtered_snv_merged_df, how="inner", left_index=True, right_index=True)

    # Heatmap after excluding mutations shared in all samples
    fig, ax = plt.subplots(figsize=(8,15))
    fig.suptitle(f"Heatmap of number {pt} patient")
    sns.heatmap(data = filtered_snv_merged_df, cbar=False, cmap=["white","steelblue"], ax=ax)
    plt.show()
    without_shared_count = filtered_snv_merged_df.shape
    print(f"The total number of merged SNVs excluding mutations shared in all samples in patient {pt} is {without_shared_count} !!!")
    
    # Saving for phylogeny tree
    snv_merged[pt] = filtered_snv_merged_df



snv_merged_forphylo = copy.deepcopy(snv_merged)

for pt in dfs_snv_filtered_3rd.keys():
    _df = snv_merged_forphylo[pt]
    _df.reset_index(drop=True, inplace=True)
    _df = _df.transpose()
    _df = _df.replace({True:1,False:0})
    snv_merged_forphylo[pt] = _df


for pt in dfs_snv_filtered_3rd.keys():
    _df = snv_merged_forphylo[pt]
    rows = _df.shape[0]
    cols = _df.shape[1]
    
    # Patient 70 has different nomenclature pattern (ouput - )
    if pt == '70':
        _df.rename(index={k: "-".join([k.split("-")[0], k.split("-")[2]]) for k in list(_df.index)}, inplace=True)
    
    elif pt == '57':
        _df.rename(index={k: "-Ov-".join([k.split("-")[0], k.split("-")[2]]) for k in list(_df.index) if "Ovary" in k}, inplace=True)
    
    else:
        pass
   
    # Write proper input file for phylip
    # Name of samples has to be equal to or less than 10
    f = open(f"02_vcf/phylogeny/{pt}/{pt}.infile", "w")
    
    f.write(str(rows+1) + " " + str(cols) + "\n")
    f.write("dummy".ljust(10) + cols*"0" + "\n")
    for sample in list(_df.index):
        row = _df.loc[sample,:]
        row_string_withenter = row.to_string(index=False)
        row_list = list(row_string_withenter)
    
        for string in row_list:
            if string == "\n":
                row_list.remove(string)
            
        f.write(str(sample).ljust(10) + " " + str("".join(row_list)) + "\n")
        
    f.close()

################################################################################3


 # Use phylip in the terminal, with options; jumble=100, outgroup=1 (use dummy)


############################################################################


# Counting mutations per branch


muts_br = defaultdict(dict)
for pt in snv_merged.keys():

    counts_df = copy.deepcopy(snv_merged[pt])
    t=ete3.Tree(f"/home/users/wonhee720/Projects/03_CRC/02_vcf/phylogeny/{pt}/outtree")
    nodes_ls = []
    for n in t.traverse():
        if pt == '70':
            pt_70_names = n.get_leaf_names()
            pt_70_names_new = []
            for name in pt_70_names:
                if name == 'dummy':
                    pt_70_names_new.append(name)
                else:
                    temp_ls = name.split("-")
                    temp_ls.append("Colon")
                    order = [0,2,1]
                    temp_ls = [temp_ls[i] for i in order]
                    pt_70_names_new.append("-".join(temp_ls))
            n.name = pt_70_names_new
            nodes_ls.append(list(n.name))
        
        
        elif pt == '57':
            pt_57_names = n.get_leaf_names()
            pt_57_names_new = []
            for name in pt_57_names:
                if 'Ov' in name:
                    name = name.replace("Ov", "Ovary")
                    pt_57_names_new.append(name)
                else:
                    pt_57_names_new.append(name)
            n.name = pt_57_names_new
            nodes_ls.append(list(n.name))
            
        else:
            n.name = n.get_leaf_names()   
            nodes_ls.append(list(n.name))
    
    for i, node in enumerate(nodes_ls):
        # Get a list containing all samples except dummy
        if i ==0:
            root = copy.deepcopy(node)
            root.remove("dummy")
            
        # pass dummy nodes
        elif "dummy" in node:
            pass
        
        else:
            # Make a copy of a disposable dataframe
            _df = copy.deepcopy(snv_merged[pt])
            yes = []
            no = []
            
            for x in root:
                # Get list containing samples that belong to this node
                if x in node:
                    yes.append(x)
                # Get list containing samples that doesnt belong to this node
                else:
                    no.append(x)
            
            # Filter the copied dataframe using "yes" and "no" lists
            for y_comp in yes:
                _df = _df[_df[y_comp]==1]
            
            for n_comp in no:
                _df = _df[_df[n_comp]!=1]
            
            # Save it to a mutation per branch dictionary
            muts_br[pt][str(node)] = _df
            
    for n in t.traverse():
        node = n.name
        n.name = ", ".join(node)
        if "dummy" in node:
            # Remove dummy nodes
            n.delete()
        else:
            # Save mutation counts to each branches
            n.dist = muts_br[pt][str(node)].shape[0]
            print(node)
            print(n.dist)
            
    # Format 5 = "internal and leaf branches + leaf names"
    t.write(format=5, outfile=f"/home/users/wonhee720/Projects/03_CRC/02_vcf/phylogeny/{pt}/outtree_withlength")
    
    # Format 3 = "all branches + all names"
    t.write(format=3, outfile=f"/home/users/wonhee720/Projects/03_CRC/02_vcf/phylogeny/{pt}/outtree_withlength_allname")
    
    # Read the new newick file (with branch length)
    tree = Phylo.read(f"/home/users/wonhee720/Projects/03_CRC/02_vcf/phylogeny/{pt}/outtree_withlength", "newick")
    
    # Draw tree using Phylo
    tree.ladderize() 
    Phylo.draw(tree)


# Get signature of the root

fasta = pysam.FastaFile("/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta")

for pt, df in snv_root.items():
    df.reset_index(inplace=True)
    df.sort_values(by=["CHROM", "POS"], inplace=True)
    df[['context_3', 'substitution']] = df.apply(lambda x: scripts.sig.context_3(fasta, x['CHROM'], x['POS'], x['ALT']), axis=1, result_type='expand')


sig_profile_br = defaultdict(dict)

for pt, df in snv_root.items():
        print("signature in " + "root of patient " + str(pt) + "!!!")
        result = scripts.sig.sample_decomposition(scripts.sig.data_generation_96(df=df, column_name='context_3', sample_id=str(pt)), cancer="Colon")
        display(result[0])
        
        # Specify the root name and save signature profile to dictionary
        root = list(muts_br[pt].keys())[0]
        sig_profile_br[pt][root] = result[2]


# Save the number of mutation in root to python file containing that variable
for pt, df in snv_root.items():
    root_count = df.shape[0]
    f = open(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs/{pt}/root.py", "w")
    f.write("NUM = ")
    f.write(str(root_count))
    f.close()



# Get signature per branch

for pt in snv_merged.keys():
    for k, v in muts_br[pt].items():
        v.reset_index(inplace=True)
        v.sort_values(by=["CHROM","POS"], inplace=True)



for pt in snv_merged.keys():
    for k,v in muts_br[pt].items():
        if v.shape[0] == 0:
            pass
        else:
            v[['context_3', 'substitution']] = v.apply(lambda x: scripts.sig.context_3(fasta, x['CHROM'], x['POS'], x['ALT']), axis=1, result_type='expand')
            muts_br[pt][k] = v



for pt in snv_merged.keys():
    for k,v in muts_br[pt].items():
        if v.shape[0] == 0:
            pass
        else:
            print("signature in " + "branch No " + str(k) + "!!!")
            result = scripts.sig.sample_decomposition(scripts.sig.data_generation_96(df=v, column_name='context_3', sample_id=str(k)), cancer='Colon')
            display(result[0])
            
            # Save the signature profile to dictionary
            sig_profile_br[pt][k] = result[2]




# Save the branch signature profile to a png file

for pt in sig_profile_br.keys():
    br_key = {}
    for i, (k,v) in enumerate(sig_profile_br[pt].items()): # The first item must be the ROOT !!! (due to tree_python.py script)
        tmp_dict = {}
        br_key[k] = i
        sigs = v["Signature"]
        for sig in sigs:
            muts = int(v[v["Signature"] ==  sig]['exposure'])
            tmp_dict[sig] = [muts]
            
        _df = pd.DataFrame(tmp_dict, index=["Signatures"])
        _df.plot(kind='barh', stacked=True, legend=False, color=sig_color_dict_modi)
        plt.axis('off')
        plt.savefig(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs/{pt}/{i}.png", format='png')
    
    f = open(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs/{pt}/key.py", "w")
    f.write("KEY = ")
    f.write(str(br_key))
    f.close()


# Save each root and branch 


# Save the branch signature profile to a png file

for pt in snv_merged.keys():
    for k, v in muts_br[pt].items():
        if v.shape[0] == 0:
            pass
        else:
            v.to_csv(f"/home/users/wonhee720/Projects/03_CRC/06_results/branch/{pt}/{k}.branch.snps.tsv", sep="\t", index=False)

for pt, df in snv_root.items():
    _df = copy.deepcopy(df)
    _df.set_index(['CHROM','POS','REF','ALT'], inplace=True)
    _df[sample_id_dict[pt]].to_csv(f"/home/users/wonhee720/Projects/03_CRC/06_results/branch/{pt}/{pt}.root.snps.tsv", sep="\t")


for pt in snv_merged.keys():
    for k, v in muts_br[pt].items():
        v.reset_index(inplace=True)
        v.sort_values(by=["CHROM","POS"], inplace=True)

for pt in snv_merged.keys():
    for k,v in muts_br[pt].items():
        if v.shape[0] == 0:
            pass
        else:
            v[['context_3', 'substitution']] = v.apply(lambda x: scripts.sig.context_3(fasta, x['CHROM'], x['POS'], x['ALT']), axis=1, result_type='expand')
            muts_br[pt][k] = v

for pt in snv_merged.keys():
    for k,v in muts_br[pt].items():
        if v.shape[0] == 0:
            pass
        else:
            print("signature in " + "branch No " + str(k) + "!!!")
            result = scripts.sig.sample_decomposition(scripts.sig.data_generation_96(df=v, column_name='context_3', sample_id=str(k)), cancer='Colon')
            display(result[0])
            
            # Save the signature profile to dictionary
            sig_profile_br[pt][k] = result[2]


# Save the branch signature profile to a png file

for pt in sig_profile_br.keys():
    br_key = {}
    for i, (k,v) in enumerate(sig_profile_br[pt].items()): # The first item must be the ROOT !!! (due to tree_python.py script)
        tmp_dict = {}
        br_key[k] = i
        sigs = v["Signature"]
        for sig in sigs:
            muts = int(v[v["Signature"] ==  sig]['exposure'])
            tmp_dict[sig] = [muts]
            
        _df = pd.DataFrame(tmp_dict, index=["Signatures"])
        _df.plot(kind='barh', stacked=True, legend=False, color=sig_color_dict)
        plt.axis('off')
        plt.savefig(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs/{pt}/{i}.png", format='png')
    
    f = open(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs/{pt}/key.py", "w")
    f.write("KEY = ")
    f.write(str(br_key))
    f.close()

# Checking validity of the phylogenetic tree 

for pt, v in snv_merged.items():
    muts_fittree_only = pd.concat(muts_br[pt].values())
    total_fit_counts = muts_fittree_only.shape[0]
    original_counts = v.shape[0]
    pct = round(total_fit_counts/original_counts, 4) * 100
    
    # Number of mutations that fit the tree
    print("Patient " + pt + " has " + str(total_fit_counts) + " mutations that fit into the tree !!!")
    print("It is " + str(pct) + "% of total merged SNV mutations !!!")
    
    # Heatmap of mutations that fit the tree
    fig2, ax2 = plt.subplots(figsize=(8,15))
    fig2.suptitle(f"Heatmap of number {pt} patient")
    
    sorted_muts_fittree_only = merged_for_index[pt][[]].merge(muts_fittree_only.set_index(
        ["CHROM","POS","REF","ALT"]), how='inner', left_index=True, right_index=True)
    
    sorted_muts_fittree_only.drop(columns=["context_3","substitution"], inplace=True)
    
    sns.heatmap(data = sorted_muts_fittree_only, cmap=["White","steelblue"], cbar=False, ax=ax2)
    plt.show()


















