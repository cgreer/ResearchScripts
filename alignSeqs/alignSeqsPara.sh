#!/bin/bash

fN=$1
dbName=$2
wordSize=$3
outDir=$4
numParts=$5
maxMismatches=$6

#create folder if it doesn't exist
echo checking if directory exists
if [ ! -d "${outDir}" ]; then
	echo ....nope, making it
	mkdir $outDir
fi

#remove existing files in out folder
echo Removing old stuff if stuff already exists
rm ${outDir}/*

#split the file up and put in folder
echo splitting file into directory
$HOME/exec/splitFile.sh $fN $outDir $numParts

#align each file seperately
echo submitting aligning jobs
for i in ${outDir}/*.splitted.*
do
	$HOME/exec/qJobX5.sh $HOME/scripts/alignSeqs/qCGAlign.sh $i $dbName $wordSize $i.aligned $maxMismatches
done

if false
then

#now send finishing script...
#$HOME/exec/qJob.sh $HOME/scripts/alignSeqs/qCGAlignFinish.sh $1 $2 $3 $4 $5 $6


#jobs finished?
echo waiting for jobs to finish...
checkExitSignal.sh $outDir $numParts

echo Gluing files together
#cat all the files into one file
for i in ${outDir}/*.aligned
do
	cat $i >> ${outDir}/all.aligned
done

#remove duplicates
echo removing duplicate alignments
python $HOME/scripts/removeDuplicates.py ${outDir}/all.aligned ${outDir}/../all.aligned.${maxMismatches}.nodups
fi
