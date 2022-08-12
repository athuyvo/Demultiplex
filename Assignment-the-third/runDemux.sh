#!/bin/env bash

#SBATCH --partition=bgmp        ### Partition (like a queue in PBS)
#SBATCH --job-name=deeemuxy    ### Job Name
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --cpus-per-task=1       ### Number of cpus needed per task
#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Node
#SBATCH --account=bgmp          ### Account used for job submission
#SBATCH --output=results-%j.out    ### File in which to store job output
#SBATCH --error=results-%j.err     ### File in which to store job error messages
#SBATCH --mail-user=anhthuyv@uoregon.edu ### Send email of job state changes

conda activate bgmp_310

r1="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
r2="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
r3="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
r4="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"
index="/projects/bgmp/anhthuyv/bioinfo/Bi622/Demultiplex/indexes.txt"
output="/projects/bgmp/anhthuyv/bioinfo/Bi622/Demultiplex/output/"

# r1="/projects/bgmp/anhthuyv/bioinfo/Bi622/Demultiplex/1000R1.fq"
# r2="/projects/bgmp/anhthuyv/bioinfo/Bi622/Demultiplex/1000R2.fq"
# r3="/projects/bgmp/anhthuyv/bioinfo/Bi622/Demultiplex/1000R3.fq"
# r4="/projects/bgmp/anhthuyv/bioinfo/Bi622/Demultiplex/1000R4.fq"

/usr/bin/time -v ./Demultiplex.py -R1 $r1 -R2 $r2 -R3 $r3 -R4 $r4 -i $index -o $output -r 363246735


