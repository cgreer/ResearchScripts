#!/bin/bash

for c in 50 60 70 80; do for m in {0..5}; do ./pipeline.runLink.sh $1/oRNA/oRNA.data.1.1.$c.$m $1/alignments/all.aligned.ids.filtered.1.1.$c.$m $1/oRNA/oRNA.data.template; done; done
