#!/bin/bash

#1 is name of file to be sorted

sort -r -s -n +1 -2 $1 | sort -r -s -n +2 -3 $1 | sort -r -n +4 -5 > $1.sorted
