#!/bin/bash

masterDir=$1
pRunDir=$2
objName=$3

if [ ! -e ${masterDir} ]
then
	mkdir ${masterDir}
else
	rm ${masterDir}/*
fi


for i in ${pRunDir}/run.*/${objName}/*; do cat $i >> ${masterDir}/`basename $i`; done
