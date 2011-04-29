#!/bin/bash

for i in 50 60; do for j in {0..3}; do python pipeline.runLinkFiltered.py runK50GU/oRNA.data.template runK50GU/all.aligned.0.5.tabs.filtered.1.1.$i.$j; cp runK50GU/oRNA.data.template runK50GU/oRNA.data.1.1.$i.$j; done; done
