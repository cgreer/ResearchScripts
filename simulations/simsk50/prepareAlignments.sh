#!/bin/bash

for i in $1/*.aligned; do echo $i; cat $i >> $1/all.aligned; done
python /home/chrisgre/scripts/endoClip/removeDuplicates.py removeDuplicates $1/all.aligned $1/all.aligned.nodup
python /home/chrisgre/scripts/alignSeqs/siRnaPredict.py truncateAlignments $1/all.aligned.nodup 0 5 $1/all.aligned.trun.nodup
