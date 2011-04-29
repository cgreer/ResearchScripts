#!/bin/bash

gawk '$16=="T"' $1 | sort -n -r +14 -15 > $1.results
python abridge.py 'abridge' $1.results > $1.results.abridged
