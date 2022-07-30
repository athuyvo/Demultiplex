# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read1 | 101bps | +33 |
| 1294_S1_L008_R2_001.fastq.gz | index1 | 8bps | +33 |
| 1294_S1_L008_R3_001.fastq.gz | index2 | 8bps | +33 |
| 1294_S1_L008_R4_001.fastq.gz | read2 | 101bps | +33 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.

    Read1:
    ![](https://github.com/athuyvo/Demultiplex/blob/11454d8e951426719f60da4a5e9fc0f55a9e2b19/Assignment-the-first/1294_S1_L008_R1_001.png)

    Read2:
    ![](https://github.com/athuyvo/Demultiplex/blob/master/Assignment-the-first/1294_S1_L008_R2_001.png)
    
    Read3:
    ![](https://github.com/athuyvo/Demultiplex/blob/master/Assignment-the-first/1294_S1_L008_R3_001.png)
    
    Read4:
    ![](https://github.com/athuyvo/Demultiplex/blob/master/Assignment-the-first/1294_S1_L008_R4_001.png)
    
    2. 101 for R1 and R4. 8bps for R2 and R4
    3. +33


What is a good quality score cutoff for index reads and biological read pairs to utilize for sample identification and downstream analysis, respectively? Justify your answer.

In most cases, a qscore of 20 is acceptable with a 1% chance of error. The quality score depends on the uniqueness and hamming distance of your indices. 

How many indexes have undetermined (N) base calls? (Utilize your command line tool knowledge. Submit the command(s) you used. CHALLENGE: use a one-line command)

```
dir="/projects/bgmp/shared/2017_sequencing/"
zcat ${dir}1294_S1_L008_R2_001.fastq.gz | sed -n 2~4p |  grep -c "N"
>3976613

zcat ${dir}1294_S1_L008_R3_001.fastq.gz | sed -n 2~4p |  grep -c "N"
3328051
```

    
## Part 2
1. Define the problem

We want to demultiplex 24 samples and find any indices that were swapped, have undetermined index-pairss and/or have indices that are below our quality score cut off.  

2. Describe output
We should get 24 output FASTQ files for matched indices for R1 and 24 FASTQ files for R2. There should also be two unknown fastq files for R1 and R2 for undetermined or low quality indices. There should also be two fastq files for R1 and R2 hopped indices. 

5. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
6. Pseudocode
7. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
