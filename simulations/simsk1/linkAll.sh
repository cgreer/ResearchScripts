#!/bin/bash

for N in {0..9}; do for i in 50 60; do for j in {0..5}; do python ~/scripts/endoClip/pipeline.simulatedLinkFiltered.py simulation.$N/oRNA.data.1.1.$i.$j simulation.$N/all.aligned.noIDs.sorted.trun.ids.filtered.$i.$j simulation.$N/oRNA.simSeqs; done; done; done
