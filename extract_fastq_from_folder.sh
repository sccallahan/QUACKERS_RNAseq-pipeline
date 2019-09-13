for folder in */
do
	echo ${folder}
	cd ${folder}
	rsync ./*.gz ./..
	# echo `ls ./*.gz`
	cd ..
done