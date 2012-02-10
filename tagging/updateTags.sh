#!/bin/bash
rootPath="/home/chrisgre/"

#TAG::find,bash::
#find all the files that are newer than the tag file
echo "Finding all recently modified files"
time find ${rootPath} -type f -size -30k -newer myTags.txt \( -name '*.py' -o -name '*.sh' -o -name '*.txt' \) > crawlFiles

#update all tags
echo "updating tags"
echo
python tagOps.py crawlTagFiles crawlFiles myTags.txt

