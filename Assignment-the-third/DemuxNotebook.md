##  DEMULTIPLEX ALGORITHM 
Packages versions:  
Python version: 3.10.5  
matplotlib: 3.5.2  
argparse: 1.4.0  
numpy: 1.23.1  
re: 2.2.1  
bioinfo: 0.3
_____
Define the problem: 

Parse through a lane of sequencing reads genenerated from the 2017 BGMP cohort. Demultiplex the samples to create 24 FASTQ files for read1 and read2 with the header containing the matched index1-index2. Find any indices that were swapped and undetermined index-pairs and after quality filtering. Quality cutoff used for script was 20.

Then need to create two FASTQ files with non-matched/unknown/low quality index pairs for read 1 and swapped indices.  

**Illumina FASTQ header explained:**


```
@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>
```


**Getting started**  
Start an interactive node:
```
srun --account=bgmp --partition=bgmp --nodes=1 --ntasks-per-node=1 --time=2:00:00 --cpus-per-task=1 --pty bash
```

FASTQ files location: 

```
/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz
/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz
/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz
```

Find how many base bps for all files using bash: 
```
zcat "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz" | head -2 | tail -1 | wc

>1       1     102

zcat "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz" | head -2 | tail -1 | wc

>1       1       9

zcat "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz" | head -2 | tail -1 | wc

>1       1       9

zcat "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz" | head -2 | tail -1 | wc
>  1       1     102
```

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read1 | 101bps | +33 |
| 1294_S1_L008_R2_001.fastq.gz | index1 | 8bps | +33 |
| 1294_S1_L008_R3_001.fastq.gz | index2 | 8bps | +33 |
| 1294_S1_L008_R4_001.fastq.gz | read2 | 101bps | +33 |

Minus 1 for read length because accounting for new line after each sequence. 
Qscore will be qscore +33 for new Illumina platforms.


Find how many sequences there are in fastq files: 

```
zcat "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz" | wc -l

>1452986940
Divide this by 4 = 363246735 sequences

zcat "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz" | wc -l
>1452986940
Divide this by 4 = 363246735 sequences
```

Create conda enviroment and install packages:
```
conda create --name bgmp_310 python=3.10
conda activate bgmp_310
conda install matplotlib
conda install numpy
```

Calculate mean score at each base position for all files files:
```
meanQscore.py -s <number of sequences> -b <basepairs> -f <filename> -o <output histogram filename>
```

Use sbatch script:

```
runQscoreR1.sh 
runQscoreR2.sh
runQscoreR3.sh
runQscoreR4.sh 
```

How many indexes have undetermined (N) base calls? 
```
dir="/projects/bgmp/shared/2017_sequencing/"
zcat ${dir}1294_S1_L008_R2_001.fastq.gz | sed -n 2~4p |  grep -c "N"
>3976613

zcat ${dir}1294_S1_L008_R3_001.fastq.gz | sed -n 2~4p |  grep -c "N"
>3328051
```

**For running test files**


```
./Demultiplex.py -R1 unittests/testR1.fq -R2 unittests/testR2.fq -R3 unittests/testR3.fq -R4 unittests/testR4.fq -i unittests/testindex.txt -o /projects/bgmp/anhthuyv/bioinfo/Bi622/Demultiplex/output/
```

See ```runDemux.sh``` for running actual files. 

Line counts of output files after demultiplexing to ensure correct number of lines that we started with:

```
wc -l *.fq

    30023496 AACAGCGA_R1.fq
    30023496 AACAGCGA_R2.fq
    27574204 ACGATCAG_R1.fq
    27574204 ACGATCAG_R2.fq
    37308936 AGAGTCCA_R1.fq
    37308936 AGAGTCCA_R2.fq
    29167556 AGGATAGC_R1.fq
    29167556 AGGATAGC_R2.fq
    33676832 ATCATGCG_R1.fq
    33676832 ATCATGCG_R2.fq
    22949060 ATCGTGGT_R1.fq
    22949060 ATCGTGGT_R2.fq
    13111352 CACTTCAC_R1.fq
    13111352 CACTTCAC_R2.fq
    19410564 CGATCGAT_R1.fq
    19410564 CGATCGAT_R2.fq
    14481868 CGGTAATC_R1.fq
    14481868 CGGTAATC_R2.fq
    60365860 CTAGCTCA_R1.fq
    60365860 CTAGCTCA_R2.fq
   116951644 CTCTGGAT_R1.fq
   116951644 CTCTGGAT_R2.fq
    22237744 GATCAAGG_R1.fq
    22237744 GATCAAGG_R2.fq
    12519952 GATCTTGC_R1.fq
    12519952 GATCTTGC_R2.fq
    22538636 GCTACTCT_R1.fq
    22538636 GCTACTCT_R2.fq
    27266940 GTAGCGTA_R1.fq
    27266940 GTAGCGTA_R2.fq
    29480956 GTCCTAAG_R1.fq
    29480956 GTCCTAAG_R2.fq
     2830960 hopped_R1.fq
     2830960 hopped_R2.fq
   242926668 TACCGGAT_R1.fq
   242926668 TACCGGAT_R2.fq
    35409160 TAGCCATG_R1.fq
    35409160 TAGCCATG_R2.fq
    36879488 TATGGCAC_R1.fq
    36879488 TATGGCAC_R2.fq
    12704556 TCGACAAG_R1.fq
    12704556 TCGACAAG_R2.fq
    37228944 TCGAGAGT_R1.fq
    37228944 TCGAGAGT_R2.fq
    14637952 TCGGATTC_R1.fq
    14637952 TCGGATTC_R2.fq
   141249416 TCTTCGAC_R1.fq
   141249416 TCTTCGAC_R2.fq
    54315024 TGTTCCGT_R1.fq
    54315024 TGTTCCGT_R2.fq
   355739172 unk_R1.fq
   355739172 unk_R2.fq
  2905973880 total
```

Results: 
unknown: 88934793, 24.483301412193008%
hopped: 707740, 0.19483726398807136%
matched: 273604202, 75.32186132381892%

Runtime: 
```	
Command being timed: "./Demultiplex.py -R1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -R2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -R3 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -R4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -i /projects/bgmp/anhthuyv/bioinfo/Bi622/Demultiplex/indexes.txt -o /projects/bgmp/anhthuyv/bioinfo/Bi622/Demultiplex/output/ -r 363246735"
User time (seconds): 3976.80
System time (seconds): 50.74
Percent of CPU this job got: 93%
Elapsed (wall clock) time (h:mm:ss or m:ss): 1:12:03
Average shared text size (kbytes): 0
Average unshared data size (kbytes): 0
Average stack size (kbytes): 0
Average total size (kbytes): 0
Maximum resident set size (kbytes): 300236
Average resident set size (kbytes): 0
Major (requiring I/O) page faults: 0
Minor (reclaiming a frame) page faults: 970259
Voluntary context switches: 49005
Involuntary context switches: 1087
Swaps: 0
File system inputs: 0
File system outputs: 0
Socket messages sent: 0
Socket messages received: 0
Signals delivered: 0
Page size (bytes): 4096
Exit status: 0
```
See ```results.txt``` for demultiplex summary.


