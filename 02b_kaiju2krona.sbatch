#!/bin/bash

#SBATCH --output=./log/%x.o%J
#SBATCH --error=./log/%x.e%J
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=8        # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem-per-cpu=6G
#SBATCH --time=10:00:00          # total run time limit (HH:MM:SS)

set -Euxo pipefail

OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

sample=$1 #"Kalmar93";
out="/home/ngp704/data/anvio_dirs/runs/"${sample}
kd="/projects/mjolnir1/data/databases/kaiju/2021-11-05/" # fungi, plasmids, viruses, nr, nr_euk, 

kaiju2krona -t ${kd}/nodes.dmp -n ${kd}/names.dmp -i ${out}"/"${sample}".kaiju_nr.out" -o ${out}"/"${sample}".kaiju_nr.out.krona"

ktImportText -o ${out}"/"${sample}.krona.html ${out}"/"${sample}".kaiju_nr.out.krona"
