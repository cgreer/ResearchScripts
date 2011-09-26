#!/bin/bash

#for i in {0..5}; do for c in 50 60 70 80; do ./pipeline.runLink.sh runK50Fix/oRNA/oRNA.data.1.1.$c.$i runK50Fix/alignmentsFix/all.aligned.ids.filtered.1.1.$c.$i runK50Fix/oRNA/oRNA.data.template; done; done
for i in {0..5}; do for c in 50 60 70 80; do ./pipeline.runLink.sh run.mm9.kidney/oRNA/oRNA.data.1.1.$c.$i run.mm9.kidney/alignments/all.aligned.ids.filtered.1.1.$c.$i run.mm9.kidney/oRNA/oRNA.data.template; done; done
