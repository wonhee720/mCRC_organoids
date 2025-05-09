#!/home/users/wonhee720/miniforge3/envs/genome_v2/bin/python

import os
import sys
import ete3
import pickle as pkl
import matplotlib.pyplot as plt
from ete3 import TreeStyle, Tree, faces
from ete3 import NodeStyle, TextFace, AttrFace, RectFace, CircleFace, ImgFace



tree_path = sys.argv[1] # Path to the tree newick file
scale_val = sys.argv[2] # get scale value (0.05~0.25 seems ok, higher values for smaller number of mutations)
patient = sys.argv[3] # Get patient number
sig_scale = int(sys.argv[4]) # get scale value for signature png image (250 seems ok, higher values for larger number of mutations)


tree = Tree(tree_path, format=3)

sys.path.append(f"/home/users/wonhee720/Projects/03_CRC/02_vcf/phylogeny/{patient}")

if "indels" in tree_path:
	sys.path.append(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs_indel_rescue/{patient}")
elif "sv" in tree_path:
	sys.path.append(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs_sv_rescue/{patient}")	
elif "L1" in tree_path:
	sys.path.append(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs_L1/{patient}")	
else:
	sys.path.append(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs_snv_rescue/{patient}")

from key import KEY
from root import NUM
from history import COLOR, SHAPE

fsize = 60 # font size
divider = 1.5 # divider to define font size for branch distance, legend

primary_names = ['31-2017-F', '42-SC-FT', '43-RT-FT', '52-R-FT','52-R2-FT','52-R3-FT', 
'57-L-FT', '57-CT-O1','57-CT-O2','57-CT-O3','57-CT-O4','57-CT-O5', '69-CT-FT', '69-CT-O1','69-CT-O3','69-CT-O4','69-CT-O5','69-CT-O7',
'69-CT-O8','69-CT-O9','69-CT-O10', '70-ColonT-FT', '70-Colon-T_O1','70-Colon-T_O2','70-Colon-T_O3','70-Colon-T_O4','70-Colon-T_O5','70-Colon-T_O6',
'70-Colon-T_O7','70-Colon-T_O8','70-Colon-T_O9','70-Colon-T_O10']

with open('/home/users/wonhee720/Projects/03_CRC/06_results/id_exchanger.pkl', 'rb') as f:
    id_exchanger = pkl.load(f)

# Tree style
ts = TreeStyle()
ts.show_leaf_name = False # Set this false in order to set the face with my own layout style 
ts.show_branch_length = False
ts.show_scale = False
ts.scale = float(scale_val)

if "indels" in tree_path: # 100 mutations per scale for indels
	if patient == '70':
		mut_unit = 5000
	else:
		mut_unit = 100

	f0 = RectFace(width=mut_unit*float(scale_val), height=10, bgcolor='black', fgcolor='white')
	f0.hz_align = 1
	ts.legend.add_face(f0, column = 0)
	ts.legend.add_face(TextFace(f"  ", ftype="Arial", fsize=fsize/divider, bold=True), column = 1)	
	ts.legend.add_face(TextFace(f"  {str(mut_unit)} mutations", ftype="Arial", fsize=fsize/divider, bold=True), column = 2)	
	ts.legend.add_face(TextFace(f"  ", ftype="Arial", fsize=fsize/divider, bold=True), column = 3)

elif "sv" in tree_path:
	f0 = RectFace(width=10*float(scale_val), height=10, bgcolor='black', fgcolor='white')
	f0.hz_align = 1
	ts.legend.add_face(f0, column = 0)
	ts.legend.add_face(TextFace(f"  ", ftype="Arial", fsize=fsize/divider, bold=True), column = 1)	
	ts.legend.add_face(TextFace(f"  10 mutations", ftype="Arial", fsize=fsize/divider, bold=True), column = 2)	
	ts.legend.add_face(TextFace(f"  ", ftype="Arial", fsize=fsize/divider, bold=True), column = 3)


elif "L1" in tree_path:
	f0 = RectFace(width=5*float(scale_val), height=10, bgcolor='black', fgcolor='white')
	f0.hz_align = 1
	ts.legend.add_face(f0, column = 0)
	ts.legend.add_face(TextFace(f"  ", ftype="Arial", fsize=fsize/divider, bold=True), column = 1)	
	ts.legend.add_face(TextFace(f"  5 mutations", ftype="Arial", fsize=fsize/divider, bold=True), column = 2)	
	ts.legend.add_face(TextFace(f"  ", ftype="Arial", fsize=fsize/divider, bold=True), column = 3)


else: # 1000 mutations per scale for snvs
	if patient == '70':
		mut_unit = 5000
	else:
		mut_unit = 1000

	f0 = RectFace(width=mut_unit*float(scale_val), height=10, bgcolor='black', fgcolor='white')
	f0.hz_align = 1
	ts.legend.add_face(f0, column = 0)
	ts.legend.add_face(TextFace(f"  ", ftype="Arial", fsize=fsize/divider, bold=True), column = 1)	
	ts.legend.add_face(TextFace(f"  {str(mut_unit)} mutations", ftype="Arial", fsize=fsize/divider, bold=True), column = 2)
	ts.legend.add_face(TextFace(f"  ", ftype="Arial", fsize=fsize/divider, bold=True), column = 3)

ts.legend.add_face(CircleFace(radius=40, color='white', style='circle'), column = 0) # for line change 
ts.legend.add_face(CircleFace(radius=40, color='white', style='circle'), column = 1) # for line change
ts.legend.add_face(CircleFace(radius=40, color='white', style='circle'), column = 2) # for line change
ts.legend.add_face(CircleFace(radius=40, color='white', style='circle'), column = 3) # for line change

f1 = CircleFace(radius=25, color='darkgray', style='circle')
f1.hz_align = 1 # align pictures to center
ts.legend.add_face(f1, column = 0)
ts.legend.add_face(TextFace(f"  ", ftype="Arial", fsize=fsize/divider, bold=True), column = 1)
ts.legend.add_face(TextFace(f"  Organoid", ftype="Arial", fsize=fsize/divider, bold=True), column = 2)
ts.legend.add_face(TextFace(f"  ", ftype="Arial", fsize=fsize/divider, bold=True), column = 3)


ts.legend.add_face(CircleFace(radius=10, color='white', style='circle'), column = 0) # for line change 
ts.legend.add_face(CircleFace(radius=10, color='white', style='circle'), column = 1) # for line change
ts.legend.add_face(CircleFace(radius=10, color='white', style='circle'), column = 2) # for line change
ts.legend.add_face(CircleFace(radius=40, color='white', style='circle'), column = 3) # for line change

f2 = RectFace(width=50, height=50, bgcolor='darkgray', fgcolor='white')
f2.hz_align = 1
ts.legend.add_face(f2, column = 0)
ts.legend.add_face(TextFace(f"  ", ftype="Arial", fsize=fsize/divider, bold=True), column = 1)
ts.legend.add_face(TextFace(f"  Frozen tissue", ftype="Arial", fsize=fsize/divider, bold=True), column = 2)
ts.legend.add_face(TextFace(f"  ", ftype="Arial", fsize=fsize/divider, bold=True), column = 3)

ts.legend.add_face(CircleFace(radius=40, color='white', style='circle'), column = 0) # for line change 
ts.legend.add_face(CircleFace(radius=40, color='white', style='circle'), column = 1) # for line change
ts.legend.add_face(CircleFace(radius=40, color='white', style='circle'), column = 2) # for line change
ts.legend.add_face(CircleFace(radius=40, color='white', style='circle'), column = 3) # for line change

f3 = CircleFace(radius=65, color='darkgray')
f3.hz_align = 1
f3_1 = RectFace(width=130, height=130, bgcolor='darkgray', fgcolor='white')
f3_1.hz_align = 1
ts.legend.add_face(f3, column = 0)
ts.legend.add_face(f3_1, column =1)
ts.legend.add_face(TextFace(f"  Primary tumor", ftype="Arial", fsize=fsize/divider, bold=True), column = 2)

ts.legend.add_face(CircleFace(radius=10, color='white', style='circle'), column = 0) # for line change 
ts.legend.add_face(CircleFace(radius=10, color='white', style='circle'), column = 1) # for line change
ts.legend.add_face(CircleFace(radius=10, color='white', style='circle'), column = 2) # for line change
ts.legend.add_face(CircleFace(radius=40, color='white', style='circle'), column = 3) # for line change

f4 = CircleFace(radius=25, color='darkgray')
f4.hz_align = 1
f4_1 = RectFace(width=50, height=50, bgcolor='darkgray', fgcolor='white')
f4_1.hz_align = 1
ts.legend.add_face(f4, column = 0)
ts.legend.add_face(f4_1, column =1)
ts.legend.add_face(TextFace(f"  Metastasis", ftype="Arial", fsize=fsize/divider, bold=True), column = 2)

ts.legend.add_face(CircleFace(radius=40, color='white', style='circle'), column = 0) # for line change 
ts.legend.add_face(CircleFace(radius=40, color='white', style='circle'), column = 1) # for line change
ts.legend.add_face(CircleFace(radius=40, color='white', style='circle'), column = 2) # for line change
ts.legend.add_face(CircleFace(radius=40, color='white', style='circle'), column = 3) # for line change

f5 = CircleFace(radius=25, color='darkblue', style='circle')
f5.hz_align	= 1
ts.legend.add_face(f5, column = 0)
ts.legend.add_face(TextFace(f"  ", ftype="Arial", fsize=fsize/divider, bold=True), column = 1)
ts.legend.add_face(TextFace(f"Samples before treatment", ftype="Arial", fsize=fsize/divider, bold=True), column = 2)

ts.legend.add_face(CircleFace(radius=10, color='white', style='circle'), column = 0) # for line change 
ts.legend.add_face(CircleFace(radius=10, color='white', style='circle'), column = 1) # for line change
ts.legend.add_face(CircleFace(radius=10, color='white', style='circle'), column = 2) # for line change
ts.legend.add_face(CircleFace(radius=40, color='white', style='circle'), column = 3) # for line change

f6 = CircleFace(radius=25, color='crimson', style='circle')
f6.hz_align = 1
ts.legend.add_face(f6, column = 0)
ts.legend.add_face(TextFace(f"  ", ftype="Arial", fsize=fsize/divider, bold=True), column = 1)
ts.legend.add_face(TextFace(f"Samples after treatment", ftype="Arial", fsize=fsize/divider, bold=True), column = 2)

ts.legend.add_face(CircleFace(radius=10, color='white', style='circle'), column = 0) # for line change 
ts.legend.add_face(CircleFace(radius=10, color='white', style='circle'), column = 1) # for line change
ts.legend.add_face(CircleFace(radius=10, color='white', style='circle'), column = 2) # for line change
ts.legend.add_face(CircleFace(radius=40, color='white', style='circle'), column = 3) # for line change


ts.legend_position = 3


# Define layout function
def my_layout(node):

	# Node style 
	ns = NodeStyle()
	ns["hz_line_width"] = 7
	ns["vt_line_width"] = 7

	if "indels" in tree_path: # For indel, some signature pngs are absent (especially for indel count around only 1)
		if node.is_root(): # For root (no face showing leaf name or branch length)
			leaf_name = node.name
			node.dist = NUM
			branch_dist = int(node.dist)

			ns['size'] = 0 # Hide root node circles

			if os.path.isfile(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs_indel_rescue/{patient}/br0.png"): # If there is no mutations in the root (patient 52), do not make face for signature
				sig_face = ImgFace(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs_indel_rescue/{patient}/br0.png", width=1.1*NUM*float(scale_val), height=120) # Need to get root number seperately because it is not included in the tree file
				dist_face = TextFace("  " + str(branch_dist) + "  \n", fsize=fsize/divider, bold=True, ftype="Arial")
				dist_face2 = TextFace("    ", fsize=fsize/divider, bold=True, ftype="Arial") 

				faces.add_face_to_node(dist_face, node, column=0, position='branch-top')
				faces.add_face_to_node(dist_face2, node, column=1, position='branch-top')			
				faces.add_face_to_node(sig_face, node, column=0, position='float')

		elif (node.is_leaf()) and ~(node.is_root()) : # For leaf nodes
			leaf_name = node.name # Leaf name
			branch_dist = int(node.dist) # Number of mutations of branches

			ns['shape'] = SHAPE[leaf_name]
			ns['fgcolor'] =  COLOR[leaf_name]

			if leaf_name in primary_names:
				ns['size'] = 130 # Size of leaves for primary samples
			else:
				ns['size'] = 50 # Size of leaves for metastasis samples

			leaf_name_list = ["'" + x + "'" for x in leaf_name.split("_ ")] # Intermediate name list
			branch_name = "[" + ", ".join(leaf_name_list) + "]" # Branch name as written in jupyter notebook (ex: "[57-CT-FT, 57-Ovary-O1]")

			
			if (KEY.get(branch_name) == None) : # Doesn't need signature image face if mutation count is 0 or br0 (root)

				name_face = TextFace(str(" " + id_exchanger[leaf_name]), fsize=fsize, bold=True, ftype="Arial") # For leaf, get a clean leaf name and show it on a textface 
				dist_face = TextFace("  " + str(branch_dist) + "  \n", fsize=fsize/divider, bold=True, ftype="Arial") # For branch, get the branch dist and show it on a textface
				dist_face2 = TextFace("    ", fsize=fsize/divider, bold=True, ftype="Arial") 

				faces.add_face_to_node(name_face, node, column=0, position='branch-right') # Add faces to node
				faces.add_face_to_node(dist_face, node, column=0, position='branch-top')
				faces.add_face_to_node(dist_face2, node, column=1, position='branch-top')

			else:
				branch_short_name = str(KEY[branch_name]) # Short name in number for each branch (too long file names cannot be written in png file)
				if (os.path.isfile(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs_indel_rescue/{patient}/{branch_short_name}.png")):
					name_face = TextFace(str(" " + id_exchanger[leaf_name]), fsize=fsize, bold=True, ftype="Arial") # For leaf, get a clean leaf name and show it on a textface 
					dist_face = TextFace("  " + str(branch_dist) + "  \n", fsize=fsize/divider, bold=True, ftype="Arial") # For branch, get the branch dist and shot it on a textface
					dist_face2 = TextFace("    ", fsize=fsize/divider, bold=True, ftype="Arial") 
					
					sig_face = ImgFace(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs_indel_rescue/{patient}/{branch_short_name}.png", width=1.1*branch_dist*float(scale_val), height=120) # Get the signature context png file per each branch
	
					faces.add_face_to_node(name_face, node, column=0, position='branch-right') # Add faces to node
					faces.add_face_to_node(dist_face, node, column=0, position='branch-top')
					faces.add_face_to_node(dist_face2, node, column=1, position='branch-top')
					faces.add_face_to_node(sig_face, node, column=0, position='float')

				else:
					name_face = TextFace(str(" " + id_exchanger[leaf_name]), fsize=fsize, bold=True, ftype="Arial") # For leaf, get a clean leaf name and show it on a textface 
					dist_face = TextFace("  " + str(branch_dist) + "  \n", fsize=fsize/divider, bold=True, ftype="Arial") # For branch, get the branch dist and show it on a textface
					dist_face2 = TextFace("    ", fsize=fsize/divider, bold=True, ftype="Arial") 

					faces.add_face_to_node(name_face, node, column=0, position='branch-right') # Add faces to node
					faces.add_face_to_node(dist_face, node, column=0, position='branch-top')
					faces.add_face_to_node(dist_face2, node, column=1, position='branch-top')


		else: # For internal nodes (no face showing leaf name)
			leaf_name = node.name
			branch_dist = int(node.dist)

			ns['size'] = 0 # Hide internal node circles
			
			leaf_name_list = ["'" + x + "'" for x in leaf_name.split("_ ")] # Intermediate name list
			branch_name = "[" + ", ".join(leaf_name_list) + "]" # Branch name as written in jupyter notebook (ex: "[57-CT-FT, 57-Ovary-O1]")

			if (KEY.get(branch_name) == None): # Doesn't need signature image face if mutation count is 0 or br0 (root)
				
				dist_face = TextFace("  " + str(branch_dist) + "  \n", fsize=fsize/divider, bold=True, ftype="Arial") 
				dist_face2 = TextFace("    ", fsize=fsize/divider, bold=True, ftype="Arial") 

				faces.add_face_to_node(dist_face, node, column=0, position='branch-top')
				faces.add_face_to_node(dist_face2, node, column=1, position='branch-top')

			
			else:
				branch_short_name = str(KEY[branch_name]) # Short name in number for each branch (too long file names cannot be written in png file)
				if (os.path.isfile(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs_indel_rescue/{patient}/{branch_short_name}.png")):
					dist_face = TextFace("  " + str(branch_dist) + "  \n", fsize=fsize/divider, bold=True, ftype="Arial") # For branch, get the branch dist and shot it on a textface
					dist_face2 = TextFace("    ", fsize=fsize/divider, bold=True, ftype="Arial") 
					
					sig_face = ImgFace(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs_indel_rescue/{patient}/{branch_short_name}.png", width=1.1*branch_dist*float(scale_val), height=120) # Get the signature context png file per each branch
	
					faces.add_face_to_node(dist_face, node, column=0, position='branch-top')
					faces.add_face_to_node(dist_face2, node, column=1, position='branch-top')
					faces.add_face_to_node(sig_face, node, column=0, position='float')

				else:
					dist_face = TextFace("  " + str(branch_dist) + "  \n", fsize=fsize/divider, bold=True, ftype="Arial") # For branch, get the branch dist and show it on a textface
					dist_face2 = TextFace("    ", fsize=fsize/divider, bold=True, ftype="Arial") 

					faces.add_face_to_node(dist_face, node, column=0, position='branch-top')
					faces.add_face_to_node(dist_face2, node, column=1, position='branch-top')
				


	elif "sv" in tree_path: # For sv, some signature pngs are absent 
		if node.is_root(): # For root (no face showing leaf name or branch length)
			leaf_name = node.name
			node.dist = NUM
			branch_dist = int(node.dist)

			ns['size'] = 0 # Hide root node circles

			if os.path.isfile(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs_sv_rescue/{patient}/br0.png"): # If there is no mutations in the root (patient 52), do not make face for signature
				sig_face = ImgFace(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs_sv_rescue/{patient}/br0.png", width=1.1*NUM*float(scale_val), height=120) # Need to get root number seperately because it is not included in the tree file
				dist_face = TextFace("  " + str(branch_dist) + "  \n", fsize=fsize/divider, bold=True, ftype="Arial")
				dist_face2 = TextFace("    ", fsize=fsize/divider, bold=True, ftype="Arial") 

				faces.add_face_to_node(dist_face, node, column=0, position='branch-top')
				faces.add_face_to_node(dist_face2, node, column=1, position='branch-top')		
				faces.add_face_to_node(sig_face, node, column=0, position='float')

		elif (node.is_leaf()) and ~(node.is_root()) : # For leaf nodes
			leaf_name = node.name # Leaf name
			branch_dist = int(node.dist) # Number of mutations of branches

			ns['shape'] = SHAPE[leaf_name]
			ns['fgcolor'] =  COLOR[leaf_name]

			if leaf_name in primary_names:
				ns['size'] = 130 # Size of leavees for primary samples
			else:
				ns['size'] = 50 # Size of leaves for metastasis samples

			leaf_name_list = ["'" + x + "'" for x in leaf_name.split("_ ")] # Intermediate name list
			branch_name = "[" + ", ".join(leaf_name_list) + "]" # Branch name as written in jupyter notebook (ex: "[57-CT-FT, 57-Ovary-O1]")
			branch_short_name = str(KEY[branch_name]) # Short name in number for each branch (too long file names cannot be written in png file)
			
			if ((not os.path.isfile(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs_sv_rescue/{patient}/{branch_short_name}.png")) | (branch_short_name == 'br0')) : # Doesn't need signature image face if mutation count is 0 or br0 (root)
				leaf_name_list = ["'" + x + "'" for x in leaf_name.split("_ ")] # Intermediate name list

				name_face = TextFace(str(" "+id_exchanger[leaf_name]), fsize=fsize, bold=True, ftype="Arial") # For leaf, get a clean leaf name and show it on a textface 
				dist_face = TextFace("  " + str(branch_dist) + "  \n", fsize=fsize/divider, bold=True, ftype="Arial") # For branch, get the branch dist and show it on a textface
				dist_face2 = TextFace("    ", fsize=fsize/divider, bold=True, ftype="Arial") 

				faces.add_face_to_node(name_face, node, column=0, position='branch-right') # Add faces to node
				faces.add_face_to_node(dist_face, node, column=0, position='branch-top')
				faces.add_face_to_node(dist_face2, node, column=1, position='branch-top')


			else:
				name_face = TextFace(str(" "+id_exchanger[leaf_name]), fsize=fsize, bold=True, ftype="Arial") # For leaf, get a clean leaf name and show it on a textface 
				dist_face = TextFace("  " + str(branch_dist) + "  \n", fsize=fsize/divider, bold=True, ftype="Arial") # For branch, get the branch dist and shot it on a textface
				dist_face2 = TextFace("    ", fsize=fsize/divider, bold=True, ftype="Arial") 
				
				sig_face = ImgFace(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs_sv_rescue/{patient}/{branch_short_name}.png", width=1.1*branch_dist*float(scale_val), height=120) # Get the signature context png file per each branch

				faces.add_face_to_node(name_face, node, column=0, position='branch-right') # Add faces to node
				faces.add_face_to_node(dist_face, node, column=0, position='branch-top')
				faces.add_face_to_node(dist_face2, node, column=1, position='branch-top')
				faces.add_face_to_node(sig_face, node, column=0, position='float')


		else: # For internal nodes (no face showing leaf name)
			leaf_name = node.name
			branch_dist = int(node.dist)

			ns['size'] = 0 # Hide internal node circles
			
			leaf_name_list = ["'" + x + "'" for x in leaf_name.split("_ ")] # Intermediate name list
			branch_name = "[" + ", ".join(leaf_name_list) + "]" # Branch name as written in jupyter notebook (ex: "[57-CT-FT, 57-Ovary-O1]")
			branch_short_name = str(KEY[branch_name]) # Short name in number for each branch (too long file names cannot be written in png file)
			
			if ((not os.path.isfile(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs_sv_rescue/{patient}/{branch_short_name}.png")) | (branch_short_name == 'br0')): # Doesn't need signature image face if mutation count is 0 or br0 (root)
				leaf_name_list = ["'" + x + "'" for x in leaf_name.split("_ ")] # Intermediate name list

				dist_face = TextFace("  " + str(branch_dist) + "  \n", fsize=fsize/divider, bold=True, ftype="Arial") 
				dist_face2 = TextFace("    ", fsize=fsize/divider, bold=True, ftype="Arial") 

				faces.add_face_to_node(dist_face, node, column=0, position='branch-top')
				faces.add_face_to_node(dist_face2, node, column=1, position='branch-top')

			
			else:
				sig_face = ImgFace(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs_sv_rescue/{patient}/{branch_short_name}.png", width=1.1*branch_dist*float(scale_val), height=120)
				dist_face = TextFace("  " + str(branch_dist) + "  \n", fsize=fsize/divider, bold=True, ftype="Arial")
				dist_face2 = TextFace("    ", fsize=fsize/divider, bold=True, ftype="Arial") 

				faces.add_face_to_node(dist_face, node, column=0, position='branch-top')
				faces.add_face_to_node(dist_face2, node, column=1, position='branch-top')
				faces.add_face_to_node(sig_face, node, column=0, position='float')
	

	elif "L1rate" in tree_path: # For L1rate, no pngs require
		if node.is_root(): # For root (no face showing leaf name or branch length)
			leaf_name = node.name
			#node.dist = NUM
			branch_dist = round(float(node.dist),1)

			ns['size'] = 0 # Hide root node circles

			if NUM!=0: # If there is no mutations in the root (patient 52), do not make face for signature
				dist_face = TextFace("  " + str(branch_dist) + "  \n", fsize=fsize/divider, bold=True, ftype="Arial")
				dist_face2 = TextFace("    ", fsize=fsize/divider, bold=True, ftype="Arial") 

				faces.add_face_to_node(dist_face, node, column=0, position='branch-top')
				faces.add_face_to_node(dist_face2, node, column=1, position='branch-top')		
				
		elif (node.is_leaf()) and ~(node.is_root()) : # For leaf nodes
			leaf_name = node.name # Leaf name
			branch_dist = round(float(node.dist),1) # Number of mutations of branches

			ns['shape'] = SHAPE[leaf_name]
			ns['fgcolor'] =  COLOR[leaf_name]

			if leaf_name in primary_names:
				ns['size'] = 130 # Size of leavees for primary samples
			else:
				ns['size'] = 50 # Size of leaves for metastasis samples

			leaf_name_list = ["'" + x + "'" for x in leaf_name.split("_ ")] # Intermediate name list
			branch_name = "[" + ", ".join(leaf_name_list) + "]" # Branch name as written in jupyter notebook (ex: "[57-CT-FT, 57-Ovary-O1]")
			branch_short_name = str(KEY[branch_name]) # Short name in number for each branch (too long file names cannot be written in png file)
			
			leaf_name_list = ["'" + x + "'" for x in leaf_name.split("_ ")] # Intermediate name list

			name_face = TextFace(str(" "+id_exchanger[leaf_name]), fsize=fsize, bold=True, ftype="Arial") # For leaf, get a clean leaf name and show it on a textface 
			dist_face = TextFace("  " + str(branch_dist) + "  \n", fsize=fsize/divider, bold=True, ftype="Arial") # For branch, get the branch dist and show it on a textface
			dist_face2 = TextFace("    ", fsize=fsize/divider, bold=True, ftype="Arial") 

			faces.add_face_to_node(name_face, node, column=0, position='branch-right') # Add faces to node
			faces.add_face_to_node(dist_face, node, column=0, position='branch-top')
			faces.add_face_to_node(dist_face2, node, column=1, position='branch-top')

		else: # For internal nodes (no face showing leaf name)
			leaf_name = node.name
			branch_dist = round(float(node.dist),1)

			ns['size'] = 0 # Hide internal node circles
			
			leaf_name_list = ["'" + x + "'" for x in leaf_name.split("_ ")] # Intermediate name list
			branch_name = "[" + ", ".join(leaf_name_list) + "]" # Branch name as written in jupyter notebook (ex: "[57-CT-FT, 57-Ovary-O1]")
			branch_short_name = str(KEY[branch_name]) # Short name in number for each branch (too long file names cannot be written in png file)
			
			leaf_name_list = ["'" + x + "'" for x in leaf_name.split("_ ")] # Intermediate name list

			dist_face = TextFace("  " + str(branch_dist) + "  \n", fsize=fsize/divider, bold=True, ftype="Arial") 
			dist_face2 = TextFace("    ", fsize=fsize/divider, bold=True, ftype="Arial") 

			faces.add_face_to_node(dist_face, node, column=0, position='branch-top')
			faces.add_face_to_node(dist_face2, node, column=1, position='branch-top')

	
	elif "L1" in tree_path: # For L1, some signature pngs are absent 
		if node.is_root(): # For root (no face showing leaf name or branch length)
			leaf_name = node.name
			node.dist = NUM
			branch_dist = int(node.dist)

			ns['size'] = 0 # Hide root node circles

			if os.path.isfile(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs_L1/{patient}/br0.png"): # If there is no mutations in the root (patient 52), do not make face for signature
				sig_face = ImgFace(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs_L1/{patient}/br0.png", width=1.1*NUM*float(scale_val), height=120) # Need to get root number seperately because it is not included in the tree file
				dist_face = TextFace("  " + str(branch_dist) + "  \n", fsize=fsize/divider, bold=True, ftype="Arial")
				dist_face2 = TextFace("    ", fsize=fsize/divider, bold=True, ftype="Arial") 

				faces.add_face_to_node(dist_face, node, column=0, position='branch-top')
				faces.add_face_to_node(dist_face2, node, column=1, position='branch-top')		
				faces.add_face_to_node(sig_face, node, column=0, position='float')

		elif (node.is_leaf()) and ~(node.is_root()) : # For leaf nodes
			leaf_name = node.name # Leaf name
			branch_dist = int(node.dist) # Number of mutations of branches

			ns['shape'] = SHAPE[leaf_name]
			ns['fgcolor'] =  COLOR[leaf_name]

			if leaf_name in primary_names:
				ns['size'] = 130 # Size of leavees for primary samples
			else:
				ns['size'] = 50 # Size of leaves for metastasis samples

			leaf_name_list = ["'" + x + "'" for x in leaf_name.split("_ ")] # Intermediate name list
			branch_name = "[" + ", ".join(leaf_name_list) + "]" # Branch name as written in jupyter notebook (ex: "[57-CT-FT, 57-Ovary-O1]")
			branch_short_name = str(KEY[branch_name]) # Short name in number for each branch (too long file names cannot be written in png file)
			
			if ((not os.path.isfile(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs_L1/{patient}/{branch_short_name}.png")) | (branch_short_name == 'br0')) : # Doesn't need signature image face if mutation count is 0 or br0 (root)
				leaf_name_list = ["'" + x + "'" for x in leaf_name.split("_ ")] # Intermediate name list

				name_face = TextFace(str(" "+id_exchanger[leaf_name]), fsize=fsize, bold=True, ftype="Arial") # For leaf, get a clean leaf name and show it on a textface 
				dist_face = TextFace("  " + str(branch_dist) + "  \n", fsize=fsize/divider, bold=True, ftype="Arial") # For branch, get the branch dist and show it on a textface
				dist_face2 = TextFace("    ", fsize=fsize/divider, bold=True, ftype="Arial") 

				faces.add_face_to_node(name_face, node, column=0, position='branch-right') # Add faces to node
				faces.add_face_to_node(dist_face, node, column=0, position='branch-top')
				faces.add_face_to_node(dist_face2, node, column=1, position='branch-top')


			else:
				name_face = TextFace(str(" "+id_exchanger[leaf_name]), fsize=fsize, bold=True, ftype="Arial") # For leaf, get a clean leaf name and show it on a textface 
				dist_face = TextFace("  " + str(branch_dist) + "  \n", fsize=fsize/divider, bold=True, ftype="Arial") # For branch, get the branch dist and shot it on a textface
				dist_face2 = TextFace("    ", fsize=fsize/divider, bold=True, ftype="Arial") 
				
				sig_face = ImgFace(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs_L1/{patient}/{branch_short_name}.png", width=1.1*branch_dist*float(scale_val), height=120) # Get the signature context png file per each branch

				faces.add_face_to_node(name_face, node, column=0, position='branch-right') # Add faces to node
				faces.add_face_to_node(dist_face, node, column=0, position='branch-top')
				faces.add_face_to_node(dist_face2, node, column=1, position='branch-top')
				faces.add_face_to_node(sig_face, node, column=0, position='float')


		else: # For internal nodes (no face showing leaf name)
			leaf_name = node.name
			branch_dist = int(node.dist)

			ns['size'] = 0 # Hide internal node circles
			
			leaf_name_list = ["'" + x + "'" for x in leaf_name.split("_ ")] # Intermediate name list
			branch_name = "[" + ", ".join(leaf_name_list) + "]" # Branch name as written in jupyter notebook (ex: "[57-CT-FT, 57-Ovary-O1]")
			branch_short_name = str(KEY[branch_name]) # Short name in number for each branch (too long file names cannot be written in png file)
			
			if ((not os.path.isfile(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs_L1/{patient}/{branch_short_name}.png")) | (branch_short_name == 'br0')): # Doesn't need signature image face if mutation count is 0 or br0 (root)
				leaf_name_list = ["'" + x + "'" for x in leaf_name.split("_ ")] # Intermediate name list

				dist_face = TextFace("  " + str(branch_dist) + "  \n", fsize=fsize/divider, bold=True, ftype="Arial") 
				dist_face2 = TextFace("    ", fsize=fsize/divider, bold=True, ftype="Arial") 

				faces.add_face_to_node(dist_face, node, column=0, position='branch-top')
				faces.add_face_to_node(dist_face2, node, column=1, position='branch-top')

			
			else:
				sig_face = ImgFace(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs_L1/{patient}/{branch_short_name}.png", width=1.1*branch_dist*float(scale_val), height=120)
				dist_face = TextFace("  " + str(branch_dist) + "  \n", fsize=fsize/divider, bold=True, ftype="Arial")
				dist_face2 = TextFace("    ", fsize=fsize/divider, bold=True, ftype="Arial") 

				faces.add_face_to_node(dist_face, node, column=0, position='branch-top')
				faces.add_face_to_node(dist_face2, node, column=1, position='branch-top')
				faces.add_face_to_node(sig_face, node, column=0, position='float')



	else:
		if node.is_root(): # For root (no face showing leaf name or branch length)
			leaf_name = node.name
			node.dist = NUM
			branch_dist = int(node.dist)

			ns['size'] = 0 # Hide root node circles

			if os.path.isfile(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs_snv_rescue/{patient}/0.png"): # If there is no mutations in the root (patient 52), do not make face for signature
				sig_face = ImgFace(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs_snv_rescue/{patient}/0.png", width=1.1*NUM*float(scale_val), height=120) # Need to get root number seperately because it is not included in the tree file
				dist_face = TextFace(("  " + str(branch_dist) + "  \n"), fsize=fsize/divider, bold=True, ftype="Arial") # For branch, get the branch dist and shot it on a textface
				dist_face2 = TextFace("    ", fsize=fsize/divider, bold=True, ftype="Arial") 

				faces.add_face_to_node(dist_face, node, column=0, position='branch-top')
				faces.add_face_to_node(dist_face2, node, column=1, position='branch-top')
				faces.add_face_to_node(sig_face, node, column=0, position='float')

		elif (node.is_leaf()) and ~(node.is_root()) : # For leaf nodes
			leaf_name = node.name # Leaf name
			branch_dist = int(node.dist) # Number of mutations of branches

			ns['shape'] = SHAPE[leaf_name]
			ns['fgcolor'] =  COLOR[leaf_name]

			if leaf_name in primary_names:
				ns['size'] = 130 # Size of leavees for primary samples
			else:
				ns['size'] = 50 # Size of leaves for metastasis samples

			if branch_dist == 0: # Doesn't need signature image face if mutation count is 0
				leaf_name_list = ["'" + x + "'" for x in leaf_name.split("_ ")] # Intermediate name list

				name_face = TextFace(str(" "+id_exchanger[leaf_name]), fsize=fsize, bold=True, ftype="Arial") # For leaf, get a clean leaf name and show it on a textface 
				dist_face = TextFace(("  " + str(branch_dist) + "  \n"), fsize=fsize/divider, bold=True, ftype="Arial") # For branch, get the branch dist and show it on a textface
				dist_face2 = TextFace("    ", fsize=fsize/divider, bold=True, ftype="Arial") 

				faces.add_face_to_node(name_face, node, column=0, position='branch-right') # Add faces to node
				faces.add_face_to_node(dist_face, node, column=0, position='branch-top')
				faces.add_face_to_node(dist_face2, node, column=1, position='branch-top')


			else:
				leaf_name_list = ["'" + x + "'" for x in leaf_name.split("_ ")] # Intermediate name list
				branch_name = "[" + ", ".join(leaf_name_list) + "]" # Branch name as a written in jupyter notebook (ex: "[57-CT-FT, 57-Ovary-O1]")
				branch_short_name = str(KEY[branch_name]) # Short name in number for each branch (too long file names cannot be written in png file)

				name_face = TextFace(str(" "+id_exchanger[leaf_name]), fsize=fsize, bold=True, ftype="Arial") # For leaf, get a clean leaf name and show it on a textface 
				dist_face = TextFace(("  " + str(branch_dist) + "  \n"), fsize=fsize/divider, bold=True, ftype="Arial") # For branch, get the branch dist and shot it on a textface
				dist_face2 = TextFace("    ", fsize=fsize/divider, bold=True, ftype="Arial") 

				sig_face = ImgFace(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs_snv_rescue/{patient}/{branch_short_name}.png", width=1.1*branch_dist*float(scale_val), height=120) # Get the signature context png file per each bracnh

				faces.add_face_to_node(name_face, node, column=0, position='branch-right') # Add faces to node
				faces.add_face_to_node(dist_face, node, column=0, position='branch-top')
				faces.add_face_to_node(dist_face2, node, column=1, position='branch-top')
				faces.add_face_to_node(sig_face, node, column=0, position='float')

		else: # For internal nodes (no face showing leaf name)
			leaf_name = node.name
			branch_dist = int(node.dist)

			ns['size'] = 0 # Hide internal node circles

			if branch_dist == 0: # Doesn't need signature image face if mutation count is 0
				leaf_name_list = ["'" + x + "'" for x in leaf_name.split("_ ")] # Intermediate name list

				dist_face = TextFace(("  " + str(branch_dist) + "  \n"), fsize=fsize/divider, bold=True, ftype="Arial") 
				dist_face2 = TextFace("    ", fsize=fsize/divider, bold=True, ftype="Arial") 

				faces.add_face_to_node(dist_face, node, column=0, position='branch-top')
				faces.add_face_to_node(dist_face2, node, column=1, position='branch-top')
			
			else:
				leaf_name_list = ["'" + x + "'" for x in leaf_name.split("_ ")]
				branch_name = "[" + ", ".join(leaf_name_list) + "]"
				branch_short_name = str(KEY[branch_name])

				sig_face = ImgFace(f"/home/users/wonhee720/Projects/03_CRC/06_results/sigs_snv_rescue/{patient}/{branch_short_name}.png", width=1.1*branch_dist*float(scale_val), height=120)
				dist_face = TextFace(("  " + str(branch_dist) + "  \n"), fsize=fsize/divider, bold=True, ftype="Arial")
				dist_face2 = TextFace("    ", fsize=fsize/divider, bold=True, ftype="Arial") 

				faces.add_face_to_node(dist_face, node, column=0, position='branch-top')
				faces.add_face_to_node(dist_face2, node, column=1, position='branch-top')
				faces.add_face_to_node(sig_face, node, column=0, position='float')

	node.set_style(ns)
	


# Apply the layout function to treestyle
ts.layout_fn = my_layout

# Save figure
outpath = tree_path + ".ete3.svg"
tree.render(outpath, tree_style=ts)
