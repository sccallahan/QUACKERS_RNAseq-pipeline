#! /bin/bash

 

set -e 

set -u

set -o pipefail

 
pwd=`pwd`
root='pwd'
index_path=$1
genome_gtf=$2
fastq_dir=$3

 

cat ${fastq_dir}/sample_sheet_STAR.txt | while read -r fq1 fq2

 

do

prefix=$(basename ${fq1} .fastq.gz)

 

JobString="

#BSUB -J ${prefix}_align_STAR_and_count

#BSUB -W 12:00

#BSUB -o ${pwd}/logs 

#BSUB -e ${pwd}/logs

#BSUB -q medium

#BSUB -n 10

#BSUB -M 32

#BSUB -R rusage[mem=32]

#BSUB -N

 

module load star/2.7.2b
module load subread/1.6.3

STAR --runThreadN 10 --sjdbGTFfile ${genome_gtf} --outFileNamePrefix ${pwd}/STAR_output/${prefix} --genomeDir ${index_path} --readFilesCommand zcat --readFilesIn ${fq1} ${fq2} --outSAMtype BAM SortedByCoordinate

featureCounts -p -T 10 -a ${genome_gtf} -o ${pwd}/featureCounts_output/${prefix}_featurecounts.txt ${pwd}/STAR_output/${prefix}Aligned.sortedByCoord.out.bam
"

echo "$JobString" > ${pwd}/${prefix}_STAR_RNA_SEQ.lsf

bsub < ${pwd}/${prefix}_STAR_RNA_SEQ.lsf

done