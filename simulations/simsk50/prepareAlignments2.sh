#!/bin/bash

rm $1/all.aligned
for i in $1/*.aligned; do echo $i; cat $i >> $1/all.aligned; done
python /home/chrisgre/scripts/endoClip/removeDuplicates.py removeDuplicates $1/all.aligned $1/all.aligned.nodup
