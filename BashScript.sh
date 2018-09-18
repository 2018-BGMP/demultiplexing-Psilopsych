#!/usr/bin/bash
#SBATCH --partition=short       ### Partition (like a queue in PBS)
#SBATCH --job-name=StatsR3      ### Job Name
#SBATCH --time=1-00:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Nod

# module load python3
module purge
module load easybuild intel/2017a Python/3.6.1; which python

/usr/bin/time python3 posStatsN.py > statsR3.txt