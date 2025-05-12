#!/bin/bash
# ------------------------------------------------------------------
# [Author] Title
#          Description
# ------------------------------------------------------------------

set -Eeuo pipefail

wd="/home/ngp704/data/anvio_dirs/relative_abundance/"
cd ${wd};

virus=$1 # Ceduovirus_2_MAGs_38.BLAST
virus_blast=${wd}"/virus_BLAST/"${virus}
# echo ${virus_blast};

# rm ${wd}"/virus_BLAST/"*_hits_blast*;

for i in `ls /projects/mjolnir1/data/databases/blastdb/20220825/nt*nsq | cut -f2 -d'.'`; do
    echo $i;
    #i=01
    arxi=${virus_blast}"_"${i}".BLAST"
    echo $arxi;
    while read a; 
    do
	# nn=$(echo ${a} | awk -F'___' '{print $1"___"$2}'); ## For  anvio HMM output
	nn=$a # For the metaflye assembly
	hits=$(grep -v 'Fields' ${arxi} | grep "${nn}" -A 2 | tail -n1 | cut -f2 -d' ');
	if [ "$hits" -eq 0 ]; 
	then
	    continue
	fi
	sed -ne "/$nn/,$ p" ${arxi} | grep -v '#' > /tmp/test_${i}_${nn};
	head -n${hits} /tmp/test_${i}_${nn} >> ${wd}"/virus_BLAST/"${nn}_hits_blast_${i}.txt;
	echo $hits;
    done < <(grep '>' ${wd}"/"${virus}".fa" | sed 's/>//g');
done;
