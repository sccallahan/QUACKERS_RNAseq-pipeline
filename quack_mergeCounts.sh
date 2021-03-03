#! /bin/bash

set -e 

set -u

set -o pipefail

 
pwd=`pwd`


#### make job ####
JobString="

#BSUB -J merge_counts

#BSUB -W 08:00

#BSUB -o ${pwd}/logs

#BSUB -e ${pwd}/logs

#BSUB -q medium

#BSUB -n 10

#BSUB -M 32

#BSUB -R rusage[mem=32]

#BSUB -N

module load parallel/20201122

# From Ming Tang https://divingintogeneticsandgenomics.rbind.io/post/merge-featurecount-table-from-rnaseq/
cd featureCounts_output
ls -1  *featurecounts.txt | parallel 'cat {} | sed '1d' | cut -f7 {} > {/.}_clean.txt' 
ls -1  *featurecounts.txt | head -1 | xargs cut -f1 > genes.txt
paste genes.txt *featurecounts_clean.txt > raw_counts_table.txt
"

echo "$JobString" > ${pwd}/merge_counts.lsf

bsub < ${pwd}/merge_counts.lsf