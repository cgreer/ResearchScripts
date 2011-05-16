#!/bin/bash

for i in C_5UTR C_EXON C_INTRON C_3UTR NC_5UTR NC_EXON NC_INTRON NC_3UTR INTER; do echo $i `grep $i $1 | wc -l`; done
