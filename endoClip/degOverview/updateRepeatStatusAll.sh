#!/bin/bash

#argument is the oRNA file

for i in 1 -1; do for j in {1..22}; do echo $i $j; python degOverview.py updateRepeatStatus $1 ../repeatMasker/repeatWigs chr$j $i; done; done
for i in 1 -1; do for j in X Y M; do echo $i $j; python degOverview.py updateRepeatStatus $1 ../repeatMasker/repeatWigs chr$j $i; done; done
