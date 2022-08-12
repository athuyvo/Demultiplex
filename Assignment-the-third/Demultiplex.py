#!/usr/bin/env python

# This script demultiplexes Illumina zipped FASTQ files with indices of 8 pairs. 
# Samples are binned into their matched indices in read1 and read2 files. 
# Any unmatched indices and indices with quality scores below 20 will be binned 
# as "unknown". Indices that are suspected to have index-hopping will be binned 
# as "hopped".

import argparse
import re 
from itertools import zip_longest
import bioinfo
import gzip
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt


# Required global variables from argparse:
# Input files for Read1, Read2, Read3, and Read4. 
# Input files for indices used to demultiplex samples. 
# Path to write output files.
# Number of FASTQ records in each file.
def get_args():
    parser = argparse.ArgumentParser(description= "input fastq files")
    parser.add_argument("-R1", "--read1", required=True, type=str, help="read1 fastq filename")
    parser.add_argument("-R2", "--read2", required=True, type=str, help="read2 fastq filename")
    parser.add_argument("-R3", "--read3", required=True, type=str, help="read3 fastq filename")
    parser.add_argument("-R4", "--read4", required=True, type=str, help="read4 fastq filename")
    parser.add_argument("-i", "--index", required=True, type=str, help="indices used filename")
    parser.add_argument("-o", "--output", required=True, type=str, help="path to write output files")
    parser.add_argument("-r", "--records", required=True, type=int, help="number of fastq records")

    return(parser.parse_args())
    
args = get_args()
R1 = args.read1
R2 = args.read2
R3 = args.read3
R4 = args.read4
index_filename = args.index
pathout = args.output
records = args.records

fastq_record = {} # dictionary of current fastq record for each fastq file with line type as keys and list of line as values
file_dict = {}  # dictionary of all write files for read 1 and read1 matched index reads 
count_dict = {} # dictionary of each record counts
swapped_dict = {}   # dictionary of swapped index-pairs counts
index_dict = {} # dictionary of indices and counts

# create a dictionary of indices used from a given index file 
# create a dictionary of write files with index names and opens them 
def create_index_dict(index_file):

    # parse through each line index file and find indices
    for line in index_file:
        index = re.findall(r'[ATGC]{8}', line.strip())

        # skip first file line and add indices to a dictionary
        # create dictionary of write files with index names and
        #   opens the files
        if (len(index) != 0):
            index_dict[index[0]] = 0

            R1file = index[0] + "_R1.fq"
            R2file = index[0] + "_R2.fq"
            fh1 = open(pathout + R1file, "w") 
            fh2 = open(pathout + R2file, "w")
            file_dict[index[0]] = [fh1, fh2]    

 # create and return a dictionary for current fastq record with 
 #  current fastq line type as key and list of fastq lines as value
def create_record_dict(line_list, line_count):
    line_type = ""
    
    if line_count % 4 == 1: 
        line_type = "header" 
    elif line_count % 4 == 2  : 
        line_type = "seq" 
        line_list[2] = reverse_complement(line_list[2])
    elif line_count % 4 == 3 :
        line_type = "plus"
    else:
        line_type = "qscore"

    fastq_record[line_type] = line_list
    return (fastq_record)

# create a dictionary with filetype as key and count of records as value
def create_count_dict():
    count_dict["unknown"] = 0
    count_dict["hopped"] = 0
    count_dict["matched"] = 0

# return reverse complement of given DNA string   
def reverse_complement(seq): 
    rev_seq = ""
    base_dict = {"A":"T","T":"A","G":"C","C":"G", "N":"N"}
    for base in seq: 
        rev_seq = base_dict[base] + rev_seq

    return (rev_seq)

# write out fastq record to given file
# header is concatenated with index1 and reverse complemented index2
# fastq record is in format of:
# header: {R1, R2, R3, R4}, sequence: {R1, R2, R3, R4}, qscore: {R1, R2, R3, R4}
def write_fastq(read1_file, read2_file):
    r1_header = fastq_record["header"][0] 
    index1 = fastq_record["seq"][1]
    index2 = fastq_record["seq"][2]
    r1_seq = fastq_record["seq"][0]
    r1_qscore = fastq_record["qscore"][0] 

    r2_header = fastq_record["header"][3] 
    r2_seq = fastq_record["seq"][3]
    r2_qscore = fastq_record["qscore"][3]

    read1_file.write(r1_header + "_" + index1 + "-" + index2 + "\n" + r1_seq + "\n+\n" + r1_qscore + "\n")  
    read2_file.write(r2_header + "_" + index1 + "-" + index2 + "\n" + r2_seq + "\n+\n" + r2_qscore + "\n")

# check quality score for given index
# return True if quality score is above 20
def check_qscore(phred):
    good_score = True

    for score in phred1:
        qscore = bioinfo.convert_phred(score)
        if (qscore < 20):
            good_score = False
            break

    return (good_score)

# calculate the percentage of records for unknown, hopped, and matched samples
# divide each count in dictionary by number of records 
# write report to an output file

def calc_report(dict, title, file):    
    out.write(title + ": \n")
    for key in dict: 
        percent = dict[key]/records * 100
        out.write(key + ": " + str(dict[key]) + " reads, " + str(round(percent, 2)) + "%\n")
    out.write("\n")


# with open(R1,"r") as R1_file, open(R2, "r") as R2_file, open(R3, "r") as R3_file, \
# open(R4, "r") as R4_file, open(index_filename, "r") as index_file, \
# open(pathout + "unk_R1.fq", "w") as unkR1_file, open(pathout + "unk_R2.fq", "w") as unkR2_file, \
# open(pathout + "hopped_R1.fq", "w") as hopR1, open(pathout + "hopped_R2.fq", "w") as hopR2:


# open and read all fastq files simultaneously
# open write files for unknown and hopped indices

with gzip.open(R1,"rt") as R1_file, gzip.open(R2, "rt") as R2_file, gzip.open(R3, "rt") as R3_file, \
gzip.open(R4, "rt") as R4_file, open(index_filename, "r") as index_file, \
open(pathout + "unk_R1.fq", "w") as unkR1_file, open(pathout + "unk_R2.fq", "w") as unkR2_file, \
open(pathout + "hopped_R1.fq", "w") as hopR1, open(pathout + "hopped_R2.fq", "w") as hopR2:

    files = [R1_file, R2_file, R3_file, R4_file] # list of fastq files 
    create_index_dict(index_file) # set of indices used
    create_count_dict() # keep track of each record 

    line_count = 0 # keep track of fastq line
    
    # compile each line from each fastq file into a list 
    for read1, read2, read3, read4 in zip_longest(*files): 
        read1, read2, read3, read4 = read1.strip(), read2.strip(), read3.strip(), read4.strip()
        reads = [read1, read2, read3, read4]
        line_count += 1

        # create a dictionary with read filenames as key and list of lines as value
        fastq_record = create_record_dict(reads, line_count)    

        # reached the end of the current fastq record
        # begin demultiplexing samples
        if line_count % 4 == 0:
            index1 = fastq_record["seq"][1]
            index2 = fastq_record["seq"][2]
           
            # check for indices not in index set
            # if yes, write out to "unknown files for read1 and read2"

            if index1 not in index_dict.keys() or index2 not in index_dict.keys():
                write_fastq(unkR1_file, unkR2_file)
                count_dict["unknown"] += 1  # update record count

            # check if index1 and reverse complement of index2 match
            elif index1 == index2: 
                   
                    # check quality score at index1 and index2
                    phred1 = fastq_record["qscore"][1]
                    phred2 = fastq_record["qscore"][2]
                    
                    # check index1 quality score
                    qscore = check_qscore(phred1)

                    # if index1 quality is good, check index2
                    if qscore == True: 
                        qscore = check_qscore(phred2)
                    
                    # write to unknown file if index quality is below 20
                    #   else write to matched index files 
                    if qscore == False:    
                        write_fastq(unkR1_file, unkR2_file)
                        count_dict["unknown"] +=1 # update record count 
                    else:
                        write_fastq(file_dict[index1][0], file_dict[index2][1])
                        count_dict["matched"] +=1 # update record count 
                        index_dict[index1] += 1 # update matched index count

            # indices are in index set, but don't match eg. hopped. 
            else:  
                write_fastq(hopR1, hopR2)
                count_dict["hopped"] +=1 # update record count

                # update swapped counts for each swapped pairs in dictionary
                if (index1 + "-" + index2) in swapped_dict:
                    swapped_dict[index1 + "-" + index2] += 1
                else:
                    swapped_dict[index1 + "-" + index2] = 1


# calculate and write out summary of results
with open(pathout + "results.txt", "w") as out: 
    calc_report(count_dict, "Results Summary", out)
    calc_report(index_dict, "Matched Dual Pairs", out)
    calc_report(swapped_dict, "Swapped Pairs", out)




