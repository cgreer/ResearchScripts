
outDir=$1
maxMismatches=$2


echo waiting for jobs to finish...
checkExitSignal.sh $outDir 100

#cat all the files into one file
for i in ${outDir}/*.aligned
do	
	echo $i
	cat $i >> ${outDir}/all.aligned
done

#remove duplicates
python removeDuplicates.py ${outDir}/all.aligned ${outDir}/../all.aligned.${maxMismatches}.nodups
