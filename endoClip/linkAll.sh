#!/bin/bash

for c in 50 60; do for m in {0..5}; do ./pipeline.runLink.sh runK50GU/oRNA.data.1.1.$c.$m runK50GU/all.aligned.tabs.filtered.1.1.$c.$m runK50GU/oRNA.data.template; done; done
