#! /bin/bash

 

set -e 

set -u

set -o pipefail

 
pwd=`pwd`
root='pwd'
index_path=$1
star_gtf=$2
fc_gtf=$3
dir=$4

 

cat sample_sheet_STAR.txt | while read -r fq1 fq2

 

do

prefix=$(basename ${fq1} .fastq.gz)

 

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

STAR --runThreadN 10 --sjdbGTFfile ${star_gtf} --outFileNamePrefix ${pwd}/STAR_output/${prefix} --genomeDir ${index_path} --readFilesCommand zcat --readFilesIn ${dir}/${fq1} ${dir}/${fq2} --outSAMtype BAM SortedByCoordinate

featureCounts -p -T 10 -a ${fc_gtf} -o ${pwd}/featureCounts_output/${prefix}_featurecounts.txt ${pwd}/STAR_output/${prefix}Aligned.sortedByCoord.out.bam
"

echo "$JobString" > ${pwd}/${prefix}_STAR_RNA_SEQ.lsf

bsub < ${pwd}/${prefix}_STAR_RNA_SEQ.lsf

done