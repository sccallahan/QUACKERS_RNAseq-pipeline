#! /bin/bash


pwd=`pwd`

cd ${pwd}/STAR_output

for bam in *.bam

do

prefix=$(basename ${bam} Aligned.sortedByCoord.out.bam)

JobString="

#BSUB -J ${prefix}_fastqc

#BSUB -W 01:00

#BSUB -o ${pwd}/logs

#BSUB -e ${pwd}/logs

#BSUB -q short

#BSUB -n 8

#BSUB -M 8192

#BSUB -R rusage[mem=8192]

#BSUB -N




# we need to source activate py351, so we have multiqc already
# module load fastqc/0.11.7 

module load fastqc/0.11.7

fastqc -t 5 -f bam -o ${pwd}/fastqc_output ${bam}
"


echo "$JobString" > ${pwd}/${prefix}_fastqc.lsf

bsub < ${pwd}/${prefix}_fastqc.lsf

done
