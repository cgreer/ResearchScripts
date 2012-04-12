#!/bin/bash

for a in 50 60 70 80; do for b in {0..5}; do for c in repeat norepeat;
do	
gawk '$27 >= 2.00' oRNA.data.1.1.$a.$b.$c.results > oRNA.data.1.1.$a.$b.$c.results.2SNR
done
done
done
