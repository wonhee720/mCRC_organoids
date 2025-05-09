#!/home/users/wonhee720/anaconda3/envs/genome_new/bin/python
import os
import sys
import ete3
import matplotlib.pyplot as plt
from ete3 import TreeStyle, Tree, faces
from ete3 import NodeStyle, TextFace, AttrFace, RectFace, CircleFace, ImgFace

t = Tree(path_to_newick, format=3)
ts = TreeStyle()
ts.mode = "c" # draw tree in circular mode
ts.show_leaf_name = True
ts.optimal_scale_level = 'full' # Actually doesn't show any effect on my trees 
ts.root_opening_factor = 0.25 # Critical factor
ts.arc_start = -180
ts.arc_span = 200
t.render("/home/users/wonhee720/Projects/03_CRC/02_vcf/phylogeny/42/outtree_withlength_rescued_allname.circular.ete3.svg", tree_style=ts)