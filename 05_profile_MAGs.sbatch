#!/bin/bash

#SBATCH --output=./log/%x.o%J
#SBATCH --error=./log/%x.e%J

#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=10        # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem-per-cpu=4G

#SBATCH --time=02:00:00          # total run time limit (HH:MM:SS)
# #SBATCH --mail-type=end          # send email when job ends
# #SBATCH --mail-user=juan.rodriguez@sund.ku.dk


OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

set -Eeuxo pipefail

cd /home/ngp704/data/anvio_dirs/relative_abundance;

MAG="/home/ngp704/data/anvio_dirs/relative_abundance/MAGs_CONTIGS.db"
sample=$1 # D02abcD12_S37
batch=$2 # D12
out="/home/ngp704/data/anvio_dirs/runs/"${batch}"/07_RELATIVE_ABUNDANCE";

anvi-profile -c ${MAG} -i ${out}/${sample}.bam --num-threads 10 -o ${out}/${sample} -W;
