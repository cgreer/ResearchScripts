#!/bin/bash

for m in {0..5}; do for c in 50 60 70 80; do gawk '$27>2.0' oRNA.data.1.1.${c}.${m}.norepeat.results > oRNA.data.1.1.${c}.${m}.norepeat.results.2SNR; done; done
