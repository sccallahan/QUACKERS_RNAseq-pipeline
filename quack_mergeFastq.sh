#! /bin/bash

 

set -e 

set -u

set -o pipefail

pwd=`pwd`

fq_dir=$1

cd ${fq_dir}

# https://stackoverflow.com/questions/39230277/having-a-job-run-only-after-all-my-previous-jobs-have-finished
# base idea taken from above link


# make array of directories
dir_list=(*/)

# create job number for bwait later
jobnum=0

for dir in ${dir_list[@]}; do
	echo ${dir}

	cd ${fq_dir}/${dir}

	# from https://www.biostars.org/p/317385/
	for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq); do


		JobString="

		#BSUB -J merge_fastq_${jobnum}_${i}

		#BSUB -W 2:00

		#BSUB -o ${pwd}/logs 

		#BSUB -e ${pwd}/logs

		#BSUB -q short

		#BSUB -n 10

		#BSUB -M 16 

		#BSUB -R rusage[mem=16]

		#BSUB -N

		cd ${fq_dir}/${dir}

		echo Merging R1
		cat "$i"_L00*_R1_001.fastq.gz > "$i"_R1.fastq.gz
		echo Merging R2
		cat "$i"_L00*_R2_001.fastq.gz > "$i"_R2.fastq.gz

		"
		echo "$JobString" > ${pwd}/merge_fastq_${jobnum}_${i}.lsf

		bsub < ${pwd}/merge_fastq_${jobnum}_${i}.lsf
	
	done

	cd ..

	# echo "$JobString" > ${pwd}/merge_fastq_${dir}.lsf

	# bsub < ${pwd}/merge_fastq_${dir}.lsf

	jobnum=$((jobnum + 1))

done

JobString2="

#BSUB -J make_samplesheet

#BSUB -W 2:00

#BSUB -o ${pwd}/logs

#BSUB -e ${pwd}/logs

#BSUB -q short

#BSUB -n 10

#BSUB -M 16

#BSUB -R rusage[mem=16]

#BSUB -N

cd ${fq_dir}

## https://stackoverflow.com/questions/35231162/bash-equivalent-for-os-walk

## find each read specifically to make sure R1 and R2 are fed into STAR correctly
find ${fq_dir} -iname '*_R1.fastq.gz' > ${fq_dir}/tmp_r1.txt
find ${fq_dir} -iname '*_R2.fastq.gz' > ${fq_dir}/tmp_r2.txt

# https://stackoverflow.com/questions/14067523/moving-every-second-row-to-a-new-column-with-awk
# have to use xargs solution -- shell expands the $0 to the script name and breaks awk solution

# xargs -n2 < ${fq_dir}/tmp_samplesheet.txt > ${fq_dir}/sample_sheet_STAR.txt

# awk '{printf "%s%s",$0,NR%2?"\t":RS}' ${fq_dir}/tmp_samplesheet.txt > ${fq_dir}/sample_sheet_STAR.txt

## just paste files together as separate columns

paste tmp_r1.txt tmp_r2.txt > sample_sheet_STAR.txt

# rm tmp_samplesheet.txt

"

echo "$JobString2" > ${pwd}/make_samplesheet.lsf

while [[ `bjobs | grep "^merge_fastq*" | wc -l` -ge 1 ]]; do
	sleep 30
done

bsub < ${pwd}/make_samplesheet.lsf
