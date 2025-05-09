import sys

fn=sys.argv[1]
inputfile=open(fn, 'r')

line=inputfile.readline()
while "#"==line[0]:
	line=inputfile.readline()
line=inputfile.readline() #skip header

ofn=fn+".covstat"
outputfile=open(ofn,"w")
outputfile.write("filename\tregion\tthroughput\tavg_depth\n")
outputfile.close()

covered_region=0
covered_depth=0
avg_depth=0
while line:
	line_split=line.rstrip().split("\t")
	this_covered_region=int(line_split[3])
	this_covered_depth=int(line_split[5])
	covered_region+=this_covered_region
	covered_depth+=this_covered_depth

	line=inputfile.readline()

avg_depth=round(covered_depth/float(covered_region),2)

with open(ofn, "a") as outputfile:
	outputfile.write(fn+"\t"+str(covered_region)+"\t"+str(covered_depth)+"\t"+str(avg_depth)+"\n")

inputfile.close()