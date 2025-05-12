#!/bin/bash

#SBATCH --output=./log/%x.o%J
#SBATCH --error=./log/%x.e%J
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=10        # cpu-cores per task (>1 if multithread tasks)
#SBATCH --mem-per-cpu=4G
#SBATCH --time=04:00:00          # total run time limit (HH:MM:SS)

set -Euxo pipefail

OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK


sample=$1 #"Kalmar93";
out="/home/ngp704/data/anvio_dirs/runs/"${sample}
kd="/projects/mjolnir1/data/databases/kaiju/2021-11-05/" # fungi, plasmids, viruses, nr, nr_euk, 
#kd="/projects/mjolnir1/data/databases/kaiju/2021_02_26_refseq" # kaiju db refseq

### If we run kaiju for anvio
anvi-get-sequences-for-gene-calls -c ${out}"/03_CONTIGS/"${sample}"-contigs.db" \
     -o ${out}"/"${sample}".gene_calls.fa"

kaiju -t ${kd}/nodes.dmp \
      -f ${kd}/nr_euk/kaiju_db_nr_euk.fmi \
      -i ${out}"/"${sample}".gene_calls.fa" \
      -o ${out}"/"${sample}".kaiju_nr.out" \
      -z 10 \
      -v

kaiju-addTaxonNames -t ${kd}/nodes.dmp \
		    -n ${kd}/names.dmp \
		    -i ${out}"/"${sample}".kaiju_nr.out" \
		    -o ${out}"/"${sample}".kaiju_nr.names" \
		    -r superkingdom,phylum,class,order,family,genus,species

anvi-import-taxonomy-for-genes -i ${out}/${sample}".kaiju_nr.names" \
    -c ${out}"/03_CONTIGS/"${sample}"-contigs.db" \
    -p kaiju \
    --just-do-it

exit;

#########
# THERE IS NO NEED ANYMORE TO USE THE HACK IN THE ANVIO WEB PAGE TO ADD TO KAIJU. #
#########
# ALL THE ABOVE WORKED
# SEE DESCRIPTION ON: https://discord.com/channels/1002537821212512296/1035835970576187432/threads/1063885802049916969
# If we want to export Kaiju taxonomy to Anvio, go with:
# **ONLY**IF USING KAIJU 1.7, WE SHOULD CHANGE THE FORMAT.
# Using the script linked from the anvio repo:
# Rscript /groups/hologenomics/juarod/Projects/metacheese/kaiju2anvio.R ${out}/${sample}".kaiju_nr.names" ${out}/${sample}".kaiju_nr-fixed.names" 
# Using the in house script
# paste <(awk -F'\t' '{print $2}' ${out}/${sample}".kaiju_nr.names") <(awk -F'\t' -v OFS='\t' '{if($4 == ""){print "NA","NA","NA","NA","NA","NA","NA"} else {print $8} }' ${out}/${sample}".kaiju_nr.names" | sed 's/;/\t/g') > ${out}/${sample}".kaiju_nr-fixed.names"
# sed -i '1i gene_callers_id\tt_domain\tt_phylum\tt_class\tt_order\tt_family\tt_genus\tt_species' ${out}/${sample}".kaiju_nr-fixed.names" &
# BACK_PID=$!
# wait $BACK_PID
