#!/bin/bash

for a in 50 60; do for b in {0..5}; do for c in repeat norepeat;
do
gawk '$10!="." && $9>1.2 && $11<6 && $12<6 && $13=="F"' oRNA.data.1.1.$a.$b.$c | sort -n -r +14 -15 > oRNA.data.1.1.$a.$b.$c.results
#python abridge.py 'abridge' oRNA.data.1.1.$a.$b.$c.results > oRNA.data.1.1.$a.$b.$c.results.abridged
done; done; done
