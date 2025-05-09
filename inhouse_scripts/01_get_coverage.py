import os
import sys
import gzip
winsize=100000 # it is fixed
thr_cov=1000

fn=sys.argv[1]
ofn=sys.argv[2]
# outdir = "/home/users/kimin/projects/03_Mucosal_Melanoma/02_vcf/depth/"
# ofn = os.path.basename(fn)
# ofn = ofn.replace(".gz","")+".100kbcov"
# ofn = outdir + ofn

outputfile=open(ofn,"w")
outputfile.write("#winsize="+str(winsize)+",thr_cov="+str(thr_cov)+"(loci with > thr_cov is not counted because it is not realistic.)\n")
outputfile.write("chr\tpos\tN\teffective\tGC\tsum_cov\taverage_cov\n")
outputfile.close()

if ".gz" in fn[-3:]:
    inputfile=gzip.open(fn, 'rb')
else:
    inputfile=open(fn, 'r')

infofn= "/home/users/jklee/Projects/00_Scripts/02_Smoothened_CNA/coverage_database/r01_human_g1k_v37.fasta.100kb"
infofile=open(infofn, 'r')
infoline=infofile.readline()#header
infoline=infofile.readline()


line=inputfile.readline()
line_split=line.rstrip().split()
this_ch=line_split[0].decode('utf-8')
this_pos=int(line_split[1])
this_period=int((this_pos-1)/winsize)
prev_ch=this_ch
prev_period=this_period

found_info=0
while infoline:
    infoline_split=infoline.rstrip().split()
    info_ch=infoline_split[0]
    info_pos=int(infoline_split[1])
    info_period=int((info_pos-1)/winsize)
    if info_ch==this_ch and info_period==prev_period:
        found_info=1
        prev_info=infoline_split
        prev_info_ch=prev_info[0]
        prev_info_pos=prev_info[1]
        break
    infoline=infofile.readline()
if found_info==0:
    print("No information was found from "+str(this_ch)+"\t"+str(this_pos))
    print("Check information file")
    infofile.close()
    input()

#chr     pos     length  N       effective       At\C    G       T       GCratio
#1       1       100000  10000   90000   26735   19924   18283   25058   0.42

sum_cov=0
nocount_region = 0 #cov >thr_cov
while line:
    line_split=line.rstrip().split()
    this_ch=line_split[0].decode('utf-8')
    this_pos=int(line_split[1])
    this_cov=int(line_split[3])
    this_period=int((this_pos-1)/winsize)
    if [this_ch, this_period] != [prev_ch, prev_period]: # new period
        #print(prev_ch+":"+str(prev_period))
        effective_len=int(infoline_split[4])
        try:
            ave_cov=round(sum_cov/float(effective_len-nocount_region),2)
        except:
            ave_cov="NA"
        with open(ofn, 'a') as outputfile:
            outputfile.write(prev_info_ch+"\t"+prev_info_pos+"\t"+prev_info[3]+"\t"+prev_info[4]+"\t"+prev_info[-1]+"\t"+str(sum_cov)+"\t"+str(ave_cov)+"\n")

        prev_ch=this_ch
        prev_period=this_period

        sum_cov=0
        ave_cov=0
        nocount_region=0
        prev_period=int((this_pos-1)/winsize)

        found_info=0
        while infoline:
            infoline_split=infoline.rstrip().split()
            info_ch=infoline_split[0]
            info_pos=int(infoline_split[1])
            info_period=int((info_pos-1)/winsize)
            if info_ch==this_ch and info_period==this_period:
                found_info=1
                prev_info=infoline_split
                prev_info_ch=prev_info[0]
                prev_info_pos=prev_info[1]
                break
            infoline=infofile.readline()
        if found_info==0:
            print("No information was found from "+this_ch+"\t"+str(this_pos))
            print("Try again...")
            infofile.close()
            infofile=open(infofn, 'r')
            infoline=infofile.readline()#header
            infoline=infofile.readline()
            while infoline:
                infoline_split=infoline.rstrip().split("\t")
                info_ch=infoline_split[0]
                info_pos=int(infoline_split[1])
                info_period=int((info_pos-1)/winsize)
                if info_ch==this_ch and info_period==this_period:
                    found_info=1
                    prev_info=infoline_split
                    prev_info_ch=prev_info[0]
                    prev_info_pos=prev_info[1]
                    break
                infoline=infofile.readline()
            if found_info==0:
                print("No information was found from "+this_ch+"\t"+str(this_pos))
                print("Check information file")
                infofile.close()
                input()

        
    if this_cov > thr_cov:
        nocount_region+=1
    else:
        if line_split[2]!="N":
            sum_cov+=this_cov
    line=inputfile.readline()

#1       10037   T       19      .........,.......,^<.   =????>>HG@?????@R=?	

effective_len=int(infoline_split[4])
try:
	ave_cov=round(sum_cov/float(effective_len-nocount_region),2)
except:
	ave_cov="NA"
with open(ofn, 'a') as outputfile:
    outputfile.write(prev_info_ch+"\t"+prev_info_pos+"\t"+prev_info[3]+"\t"+prev_info[4]+"\t"+prev_info[-1]+"\t"+str(sum_cov)+"\t"+str(ave_cov)+"\n")

inputfile.close()
infofile.close()

print("done")
