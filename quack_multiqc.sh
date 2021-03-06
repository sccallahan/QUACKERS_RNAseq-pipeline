#! /bin/bash


pwd=`pwd`

JobString="

#BSUB -J MultiQC_RNA

#BSUB -W 08:00

#BSUB -o ${pwd}/logs

#BSUB -e ${pwd}/logs

#BSUB -q medium

#BSUB -n 10

#BSUB -M 16

#BSUB -R rusage[mem=16]

#BSUB -N

module load multiqc/1.8


# we need to source activate py351, so we have multiqc already
# module load fastqc/0.11.7 
# maybe add some zcat ${file} | fastqc -f fastq -o ${pwd}/fastqc_output

cd ${pwd}/multiqc_output

# include fastqc and star qc outputs
multiqc ${pwd}/fastqc_output/ ${pwd}/STAR_output/
"


echo "$JobString" > ${pwd}/MultiQC_RNA.lsf

bsub < ${pwd}/MultiQC_RNA.lsf
