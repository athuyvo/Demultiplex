#!/usr/bin/env python

import numpy as np
import bioinfo 
import argparse
from matplotlib import pyplot as plt
import gzip

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--filename", required=True, type=str, help= "input filename")
    parser.add_argument("-o", "--output", required=True, type=str, help= "out filename")
    parser.add_argument("-b", "--basepairs", required=True, type=int, help= "read length")
    parser.add_argument("-s", "--numseqs", required=True, type=int, help= "number of sequences")
    return(parser.parse_args())

args = get_args()
file = args.filename
output = args.output
bps = args.basepairs
numseqs = args.numseqs

qlist = np.zeros(bps, dtype=np.int64) #how many bytes of memory used to store numbers
# qlist = np.zeros(bps, dtype=float)

# read zipped file as text instead of binary
with gzip.open(file,"rt") as fq: 
# with open(file,"r") as fq: 
    line_count = 0
    while True:
        
        line = fq.readline().strip()
        
        if line =="":
            break
        
        if line_count % 4 == 3: 
               
            # add quality scores at each base for each fastq record
            for i in range(0, len(line)):
                value = bioinfo.convert_phred(line[i])
                qlist[i] +=value
                
        line_count+= 1 


# mean = np.zeros(bps, dtype=float)

# for index, num in enumerate(qlist):
#     qlist[index] = num/numseqs 
#     # for index, inner in enumerate(mean):
    #     print(str(index) + "\t" + str(mean[index]))

mean = [index/numseqs for index in qlist]

x = range(0,bps)
plt.bar(x, mean)
plt.xlabel("read position")
plt.ylabel("Qscore")
plt.savefig(output)