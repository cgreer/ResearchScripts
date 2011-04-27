#!/bin/bash

sleep 5h 

for i in {0..9}; do python ~/scripts/alignSeqs/pipeline.k50SimulationFilterFlat.py simulation.$i/oRNA.data simulation.$i/all.aligned.ids; sleep 10m; done
