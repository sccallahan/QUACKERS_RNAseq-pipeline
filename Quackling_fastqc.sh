#! /bin/bash


pwd=`pwd`
dir=$1 # where your raw fastqs are located
file_type=$2 # suffix for files after the "." (e.g. ".fq.gz" would be just "fq.gz")

JobString="

#BSUB -J FastQC_MultiQC_RNA

#BSUB -W 03:00

#BSUB -o ${pwd}/logs

#BSUB -e ${pwd}/logs

#BSUB -q short

#BSUB -n 10

#BSUB -M 16384

#BSUB -R rusage[mem=16384]

#BSUB -N




# we need to source activate py351, so we have multiqc already
# module load fastqc/0.11.7 
# maybe add some zcat ${file} | fastqc -f fastq -o ${pwd}/fastqc_output

module load fastqc/0.11.7

cd ${dir}
fastqc -t 10 -f fastq -o ${pwd}/fastqc_output *.${file_type}

# mkdir multiqc_output
cd ${pwd}/multiqc_output
multiqc ${pwd}/fastqc_output/*.zip
"


echo "$JobString" > ${pwd}/FastQC_MultiQC_RNA.lsf

bsub < ${pwd}/FastQC_MultiQC_RNA.lsf
