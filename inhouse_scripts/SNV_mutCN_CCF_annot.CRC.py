#2021.8.11 modified by LWH
	# changed input arguments
	# added exception for var_read + ref_read = 0 due to WES project
	# added exception for case where var_read or ref_read is empty string ('') 

#2021.9.12 - LWH
	# added sampleNameBam function to get RG name from sample id -> need to specify the bam directory
	# Need to change sampleNameBam function input per project
	# Get sampleid as third argument
	# Ignore places where total reads are more than 2500 (consumes too much time to run this script) 

#arg1: input (tsv)
#arg2: sequenza confints_CP path


#Original script: SV_mutCN_CCF_annot.py
#2018-06-14 modification format for SNV
#2018-06-28 chrY error resolve
#2018-08-09 add exception for 'NA' d/t pcawg call
#2019-04-11 resolve logarithm error which occur at purity 1.0
#2021-02-23 column number should be provided as arguments, make script python3 compatible. change output name(.scF -> .ccf)
import pandas as pd
import math
import pysam
def nCr(n,r):
	f=math.factorial
	return math.log(f(n),10)-math.log(f(r),10)-math.log((f(n-r)),10)
#/  pow(2*math.pi*r,0.5)*pow(r/math.exp(1),r) / pow(2*math.pi*n-r,0.5)*pow(n-r/math.exp(1),n-r)

def calc_scF(tcf, wt_count, var_count, CN, mCN):
	tcf=float(tcf); wt_count=int(wt_count); var_count=int(var_count); CN=float(CN); mCN=float(mCN)
	read_depth=wt_count + var_count
	Rfraction_t=CN*tcf/(CN*tcf+(1-tcf)*2)
	vaf=var_count/float(wt_count+var_count)
	global mutCN, max_mutCN, scF
	if Rfraction_t<=0:
		mutCN=0;max_mutCN=0;scF='.'
	else:
		mutCN=vaf*CN/Rfraction_t
		if CN==2 and mCN==1:
			max_mutCN=1
			scF=mutCN
		else:
			max_binomP=0
			max_mutCN=1
			for this_CN in range(1,max(int(round(CN-mCN,0))+1,1)):
				this_binomP=0
				for readnum in range(var_count,read_depth):
					if Rfraction_t == 1:
						Rfraction_t = 0.999
					this_P1=nCr(read_depth,readnum)+readnum*math.log(Rfraction_t,10)+(read_depth-readnum)*math.log(1-Rfraction_t,10)
					this_P2=nCr(readnum,var_count)+math.log(this_CN/float(CN),10)*var_count+math.log(max(1-this_CN/float(CN),0.0001),10)*(readnum-var_count)
					this_binomP+=pow(10,this_P1+this_P2)
	
				if this_binomP>= max_binomP:
					max_binomP=this_binomP
					max_mutCN=this_CN
			scF=round(mutCN/float(max_mutCN),3)

def sampleNameBam(bamFile):
    """get @RG SM: information as sample name from BAM header"""
    bam = pysam.AlignmentFile(bamFile)
    name = bam.header['RG'][0]['SM']
    return name

import sys
fn=sys.argv[1]
print(fn)
inputfile=open(fn)
line=inputfile.readline()

ofn=fn+'-ccf' 
outputfile=open(ofn,"w")

seq_path = sys.argv[2]
name = sys.argv[3]

if name == '42-PS-FT' or name == '57-L-FT' or name == '69-CT-FT' or name == '69-LN-FT' or name == '69-LS4-FT' :
	df = pd.read_csv(fn, sep="\t", header=0)
	vaf_med = df[f"vaf_all:::{name}"].median()
	tcf = vaf_med * 2 

else:
	tcf=float(pd.read_csv(seq_path, sep="\t").iloc[0,0])

tumor = sampleNameBam(f"/home/users/wonhee720/Projects/03_CRC/01_bam_new/{name}.splitmark.realigned.recal.bam") # Change directory per project


while line:
	if line[0:6]=='#CHROM':
		line_header = line.rstrip().split("\t")
		c_ref = line_header.index(f"ref_read_all:::{tumor}") if f"ref_read_all:::{tumor}" in line_header else line_header.index(f"ref_read_all:::{tumor}") # readinfo version에 따라 _all 붙이는 여부 달라짐
		c_var = line_header.index(f"var_read_all:::{tumor}") if f"var_read_all:::{tumor}" in line_header else line_header.index(f"var_read_all:::{tumor}") # readinfo version에 따라 _all 붙이는 여부 달라짐
		c_tcn = line_header.index("CNt")
		c_mcn = line_header.index("B")

		outputfile.write(line.rstrip()+"\tmutCN\tmutCH\tCCF\n")
		line=inputfile.readline()
		continue
	
	else:
		line_split=line.rstrip().split("\t")
		if line_split[c_ref] =='' or line_split[c_var] == '':
			info_list=['.','.','.']

		else:
			wt_count1=int(float(line_split[c_ref]))
			var_count1=int(float(line_split[c_var]))
		
		if line_split[c_tcn]=='.' or line_split[c_tcn]=='NA' or (wt_count1+var_count1)>2500 or line_split[c_mcn]=='NA' or (wt_count1+var_count1)==0 :
			info_list=['.','.','.']
		else:
			seqztCN1=int(line_split[c_tcn])
			seqzmCN1=int(line_split[c_mcn])

			calc_scF(tcf, wt_count1, var_count1, seqztCN1, seqzmCN1)
			mutCN1=mutCN; max_mutCN1=max_mutCN; scF1=scF
			info_list=[str(mutCN1),str(max_mutCN1), str(scF1)]
		outputfile.write(line.rstrip()+'\t'+'\t'.join(info_list)+'\n')
	line=inputfile.readline()
print("done")

