#!/bin/bash

echo ...catting
for i in $1/*.nodups.trun; do cat $i >> $1/all.aligned; done

echo ...adding IDs
python $HOME/myLibs/blankIDs.py addIDs $1/all.aligned $1/all.aligned.ids
