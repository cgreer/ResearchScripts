#!/bin/bash

sleep 3h

for c in 50 60; do for i in {0..9}; do python ~/scripts/endoClip/pipeline.filterAlignments.py simulation.$i/all.aligned.ids simulation.$i/all.aligned.ids.filtered.1.1.$c $c; done; done
