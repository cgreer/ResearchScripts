#!/bin/bash


oFN=$1
aFN=$2
originalFN=$3

cp $originalFN $oFN
python /home/chrisgre/scripts/endoClip/pipeline.runLinkFiltered.py $1 $2
