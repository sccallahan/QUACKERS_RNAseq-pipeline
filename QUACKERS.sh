#! /bin/bash
#######################################################
## Author: Carson Callahan
## Purpose: Parent script for QUACKER
## Date: 2019-06-27
## Notes: runs child scripts for pipeline
#######################################################
dir=$1
file_type=$2
hg19_path=$3
annot_path=$4

### list child scripts
### all scripts need to be in same working directory
### `source activate py351` before running, else multiqc won't work
### Syntax: bash QUACKERS.sh [fastq directory] [file extension] [hg19 path] [gene annotation file]


# make log/error directories for job submissions
mkdir -p logs
mkdir -p fastqc_output
mkdir -p multiqc_output
mkdir -p STAR_output
mkdir -p featureCounts_output

# This script runs FastQC and MultiQC
# file extension is everything after the first ".", so ".fq" would be "fq"
module load fastqc/0.11.7
bash Quackling_fastqc.sh ${dir} ${file_type}

# This script runs STAR and featureCounts
# load modules here since we'll be submitting many jobs
module load star/2.6.0c
module load subread/1.6.3
bash Quackling_STAR.sh ${hg19_path} ${annot_path} ${dir}