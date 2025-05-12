#!/bin/bash

#SBATCH --output=./log/%x.o%J
#SBATCH --error=./log/%x.e%J

#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=32        # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem-per-cpu=8G

#SBATCH --time=90:00:00          # total run time limit (HH:MM:SS)

# #SBATCH --mail-type=end          # send email when job ends
# #SBATCH --mail-user=juan.rodriguez@sund.ku.dk

OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
cnf=$1 # config.json file
anvi-run-workflow -w metagenomics -c ${cnf}
