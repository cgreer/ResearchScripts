#!/bin/bash

for c in 50 60 70 80; do for m in {0..5}; do for i in {0..49}; do python ~/scripts/endoClip/updateSignificantSequences.py updateSignificantSequences simulation.$i/oRNA.data.1.1.$c.$m simulation.$i/all.aligned.ids.filtered.1.1.$c.$m; done; done; done
