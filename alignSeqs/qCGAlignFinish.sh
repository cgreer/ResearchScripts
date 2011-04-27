#!/bin/bash

#jobs finished?
echo waiting for jobs to finish...
#checkExitSignal.sh $4 $5

echo Gluing files together
#cat all the files into one file
for i in $4/*.aligned
do
	cat $i >> $4/all.aligned
done

#remove duplicates
echo removing duplicate alignments
python $HOME/scripts/alignSeqs/removeDuplicates.py $4/all.aligned $4/../all.aligned.$6.nodups
