#!/bin/bash

python /home/chrisgre/scripts/endoClip/removeDuplicates.py removeDuplicates $1 $1.nodups
python /home/chrisgre/scripts/endoClip/siRnaPredictFlat.py truncateAlignments $1.nodups 0 5 $1.nodups.trun
