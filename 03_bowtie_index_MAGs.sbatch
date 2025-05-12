#!/bin/bash

#SBATCH --output=./log/%x.o%J
#SBATCH --error=./log/%x.e%J
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=6       # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem-per-cpu=6G
#SBATCH --time=02:00:00          # total run time limit (HH:MM:SS)

set -Euxo pipefail

OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd /projects/mjolnir1/people/ngp704/anvio_dirs/relative_abundance/ ;

bowtie2-build MAGs.fa MAGs
