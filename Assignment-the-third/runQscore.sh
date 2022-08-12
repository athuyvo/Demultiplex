#!/bin/env bash

#SBATCH --partition=bgmp        ### Partition (like a queue in PBS)
#SBATCH --job-name=QscoreR1      ### Job Name
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --cpus-per-task=1       ### Number of cpus needed per task
#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Node
#SBATCH --account=bgmp          ### Account used for job submission
#SBATCH --output=results-%j.out    ### File in which to store job output
#SBATCH --error=results-%j.err     ### File in which to store job error messages

# conda activate bgmp_310
# conda install matplotlib
# conda install numpy

dir="/projects/bgmp/shared/2017_sequencing/"
# dir="/projects/bgmp/anhthuyv/bioinfo/Bi622/Demultiplex"

/usr/bin/time -v \
./meanQscore.py -s 363246735  -b 101 -f ${dir}1294_S1_L008_R1_001.fastq.gz -o "1294_S1_L008_R1_001.png"

# ./meanQscore.py -s 100 -b 101 -f $dir/test.fastq -o "test1.png"