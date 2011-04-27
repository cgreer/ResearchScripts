#!/bin/bash

if [ ! -n "$2" ]
then
  echo "Did not pass basename"
  exit 65
fi  

F_DIR=$1
F_BASE=$2
F_PACK=$3
F_GENOME=$4
F_MER=$5

python collectMatureFrames.py -d "${F_DIR}" -b "${F_BASE}" -n "${F_PACK}" -g "${F_GENOME}" -m "${F_MER}"
