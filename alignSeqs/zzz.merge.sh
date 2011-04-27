#!/bin/bash

masterDir=$1
pRunDir=$2
objName=$3

rm -r $masterDir
mkdir $masterDir
for i in ${pRunDir}/run.*/${objName}/*; do cat $i >> ${masterDir}/`basename $i`; done
