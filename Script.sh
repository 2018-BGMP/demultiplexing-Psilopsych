#!/usr/bin/bash
#SBATCH --partition=short       ### Partition (like a queue in PBS)
#SBATCH --job-name=DeMulti      ### Job Name
#SBATCH --time=1-00:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=28     ### Number of tasks to be launched per Nod

# module load python3
module purge
module load easybuild intel/2017a Python/3.6.1; which python

/usr/bin/time python3 DeMultiAlgoArg.py -r1 "../../../../../shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz" -i1 "../../../../../shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"\
 -r2 "../../../../../shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz" -i2 "../../../../../shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz" -t 30 > out.txt
 
 find . -name '*.fastq' | xargs wc -l