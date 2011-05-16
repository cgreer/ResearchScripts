#!/bin/bash

for t in {0..9}; do for i in 50 60; do for j in {0..5}; do python pipeline.runLinkFiltered.py /home/chrisgre/scripts/simulations/simsk50FilteredMasked/simulation.$t/oRNA.data /home/chrisgre/scripts/simulations/simsk50FilteredMasked/simulation.$t/all.aligned.ids.filtered.1.1.$i.$j; cp /home/chrisgre/scripts/simulations/simsk50FilteredMasked/simulation.$t/oRNA.data /home/chrisgre/scripts/simulations/simsk50FilteredMasked/simulation.$t/oRNA.data.1.1.$i.$j; done; done; done
