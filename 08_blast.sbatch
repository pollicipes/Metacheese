#!/bin/bash

#SBATCH --output=./log/%x.o%J
#SBATCH --error=./log/%x.e%J

#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=2        # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem-per-cpu=4G

#SBATCH --time=00:20:00          # total run time limit (HH:MM:SS)
# #SBATCH --mail-type=end          # send email when job ends
# #SBATCH --mail-user=juan.rodriguez@sund.ku.dk


OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

set -Eeuxo pipefail;

num=$1
virus=$2

wd="/home/ngp704/data/anvio_dirs/relative_abundance/"
file_name=${wd}"/"${virus}".fa" # Ceduovirus_2_MAGs.fa

virus_blast=${wd}"/virus_BLAST"

mkdir -p ${virus_blast};

blastn \
    -db /projects/mjolnir1/data/databases/blastdb/20220825/nt.${num} \
    -query ${file_name} \
    -out ${virus_blast}/${virus}"_"${num}".BLAST" \
    -max_target_seqs 2 \
    -num_threads 2 \
    -outfmt "7 stitle sseqid bitscore pident qcovs qlen slen evalue qstart qend sstart send"
