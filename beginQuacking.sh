#! /bin/bash
#######################################################
## Author: Carson Callahan
## Purpose: Parent script for QUACKERS
## Date: 2019-06-27
## Notes: runs child scripts for pipeline
#######################################################

pwd=`pwd`
me=`basename "$0"`

dir=$1
index_path=$2
star_gtf=$3
fc_gtf=$4
sleep=$5

# all scripts need to be in same working directory
# `source activate py351` before running, else multiqc won't work
# Syntax: bash beginQuacking.sh [fastq directory] [index path] [gtf file] [sleep time]

#### Some rough argument checks ####
if [ $# -eq 0 ]; then
  echo "Usage: bash ${me} [fastq directory] [index path] [star gtf file] [fc gtf file] [sleep time]"
  exit 1
elif [ $# -lt 5 ]; then
  echo "Not enough arguments!"
  echo "Usage: bash ${me} [fastq directory] [index path] [star gtf file] [fc gtf file] [sleep time]"
  exit 1
elif [ $# -gt 5 ]; then
  echo "Too many arguments!"
  echo "Usage: bash ${me} [fastq directory] [index path] [star gtf file] [fc gtf file] [sleep time]"
  exit 1
fi

#### make log/error directories for job submissions ####
mkdir -p logs
mkdir -p fastqc_output
mkdir -p multiqc_output
mkdir -p STAR_output
mkdir -p featureCounts_output


#### get number of samples ####
num=$(wc -l sample_sheet_STAR.txt | awk '{print $1}')
echo "There are ${num} samples"


#### Run STAR and featurecounts ####
# load modules here since we'll be submitting many jobs
echo "submitting STAR alignment jobs..."
module load star/2.6.0c
module load subread/1.6.3
bash quack_STAR.sh ${index_path} ${star_gtf} ${fc_gtf} ${dir}
echo "STAR jobs sumbmitted!"


##### Run FastQC ####
# run on STAR aligned bams
# file_type is everything after the first "."
# run before featurecounts since this just submits jobs

# check STAR outputs
cd ${pwd}/STAR_output
# echo "`pwd`"

# wait until the first file exists
COUNTER=1
# fix this later
# maybe shopt -s nullglob; don't forget to unset
# shopt -u nullglob after
until [[ `ls *.final.out 2>/dev/null | wc -l` -gt 0 ]]
do
	echo "waiting for first STAR output... loop ${COUNTER}"
	sleep ${sleep}
	let COUNTER=${COUNTER}+1
done

echo "first file found, waiting for additional files..."

# wait until we have 1 bam for each sample
COUNTER=1
until [[ `ls *.final.out | wc -l` -eq ${num} ]]
do
	echo "waiting for ${num} total STAR outputs... loop ${COUNTER}"
	sleep ${sleep}
	let COUNTER=${COUNTER}+1
done

# run fastqc
echo "submitting fastqc jobs..."
cd ${pwd}
module load fastqc/0.11.7
bash quack_fastqc.sh
echo "fastqc jobs sumbmitted!"


#### Merge featurecounts outputs ####
# We want to wait until all the files have had counts generated
cd ${pwd}/featureCounts_output

# wait until the first file has counts
until [[ `ls *featurecounts.txt.summary 2>/dev/null | wc -l` -gt 0 ]]
do
	echo "waiting for first featurecounts output..."
	sleep ${sleep}
done

# now we wait until the number of counts files = number of samples
until [[ `ls *featurecounts.txt.summary | wc -l` -eq ${num} ]]
do
	echo "waiting for ${num} total featurecounts outputs..."
	sleep ${sleep}
done

# submit the code to merge counts into a single table
cd ${pwd}
bash quack_mergeCounts.sh
echo "mergeCounts job submitted!"


#### Run MultiQC ####
cd ${pwd}/fastqc_output

# wait for first file
COUNTER=1
until [[ `ls *.html 2>/dev/null | wc -l` -gt 0 ]]
do
	echo "waiting for first fastqc output..."
	sleep ${sleep}
	let COUNTER=${COUNTER}+1
done

# wait until we have 1 html report for each bam/sample
COUNTER=1
until [[ `ls *.html | wc -l` -eq ${num} ]]
do
	echo "waiting for ${num} total fastqc outputs... loop ${COUNTER}"
	sleep ${sleep}
	let COUNTER=${COUNTER}+1
done

# run multiqc
cd ${pwd}
# include fastqc and star qc outputs
bash quack_multiqc.sh
echo "multiqc job submitted!"


#### clean up lsf files ####
cd ${pwd}

# make lsf log directory
cd logs
mkdir -p lsf_logs

# our last job is the count merge
# wait until a file from that job exists to move lsf files
cd ${pwd}/multiqc_output

COUNTER=1
until [[ `ls multiqc_report.html 2>/dev/null | wc -l` -eq 1 ]]
do
	echo "waiting on multiqc... loop ${COUNTER}"
	sleep ${sleep}
	let COUNTER=${COUNTER}+1
done

echo "cleaning up lsf files..."
cd ${pwd}
mv ./*.lsf logs/lsf_logs/

echo "Pipeline complete!"