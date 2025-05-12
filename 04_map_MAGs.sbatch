#!/bin/bash

#SBATCH --output=./log/%x.o%J
#SBATCH --error=./log/%x.e%J
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=8        # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem-per-cpu=6G
#SBATCH --time=08:00:00          # total run time limit (HH:MM:SS)
# #SBATCH --mail-type=end          # send email when job ends
# #SBATCH --mail-user=juan.rodriguez@sund.ku.dk

OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

set -Eeuxo pipefail

cd /home/ngp704/data/anvio_dirs/relative_abundance;

sample=$1 # D02abcD12_S37
batch=$2 # D12
fastq1=$3
fastq2=$4

MAG="MAGs"
out="/home/ngp704/data/anvio_dirs/runs/"${batch}"/07_RELATIVE_ABUNDANCE";
mkdir -p ${out};

echo -e '\n' ${sample} '\n';

bowtie2 --threads 8 \
        -x ${MAG} \
        -1 ${fastq1} \
        -2 ${fastq2} \
        --no-unal \
        -S ${out}/${sample}.sam 

# covert the resulting SAM file to a BAM file:
samtools view -F 4 -bS ${out}/${sample}.sam > ${out}/${sample}-RAW.bam;
# sort and index the BAM file:
samtools sort ${out}/${sample}-RAW.bam -o ${out}/${sample}.bam;
samtools index ${out}/${sample}.bam;
# remove temporary files:
rm ${out}/${sample}.sam ${out}/${sample}-RAW.bam;

### HOW TO RUN:
# sbatch -J map ./03_map_MAGs.sbatchh $sample $batch $fastq1 $fastq2;
