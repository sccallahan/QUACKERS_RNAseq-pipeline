#! /bin/bash

 

set -e 

set -u

set -o pipefail

 
pwd=`pwd`
root='pwd'
hg19_path=$1
annot_path=$2
dir=$3

 

cat sample_sheet_STAR.txt | while read -r fq1 fq2

 

do

prefix=$(basename ${fq1} _1.fq.gz)

 

JobString="

#BSUB -J ${prefix}_align_STAR_and_count

#BSUB -W 12:00

#BSUB -o ${pwd}/logs 

#BSUB -e ${pwd}/logs

#BSUB -q medium

#BSUB -n 10

#BSUB -M 32768 

#BSUB -R rusage[mem=32768]

#BSUB -N

 

module load star/2.6.0c
module load subread/1.6.3

STAR --runThreadN 10 --sjdbGTFfile ${hg19_path}/annot/*.gtf --outFileNamePrefix ${pwd}/STAR_output/${prefix} --genomeDir ${hg19_path}/starindex --readFilesCommand zcat --readFilesIn ${dir}/${fq1} ${dir}/${fq2} --outSAMtype BAM SortedByCoordinate

featureCounts -p -T 10 -a ${annot_path} -o ${pwd}/featureCounts_output/${prefix}_featurecounts.txt ${pwd}/STAR_output/${prefix}Aligned.sortedByCoord.out.bam

# From Ming Tang https://divingintogeneticsandgenomics.rbind.io/post/merge-featurecount-table-from-rnaseq/
cd featureCounts_output
ls -1  *featurecounts.txt | parallel 'cat {} | sed '1d' | cut -f7 {} > {/.}_clean.txt' 
ls -1  *featurecounts.txt | head -1 | xargs cut -f1 > genes.txt
paste genes.txt *featurecounts_clean.txt > raw_counts_table.txt
"

echo "$JobString" > ${pwd}/${prefix}_STAR_RNA_SEQ.lsf

bsub < ${pwd}/${prefix}_STAR_RNA_SEQ.lsf

done