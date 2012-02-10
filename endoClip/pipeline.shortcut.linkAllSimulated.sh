#!/bin/bash

#for i in {0..5}; do for c in 50 60 70 80; do ./pipeline.runLink.sh runK50Fix/oRNA/oRNA.data.1.1.$c.$i runK50Fix/alignmentsFix/all.aligned.ids.filtered.1.1.$c.$i runK50Fix/oRNA/oRNA.data.template; done; done
for s in {0..99}
do
	for i in {0..5}
	do 
		for c in 50 60 70 80 
		do
			echo sim $s $i $c 
			python pipeline.simulatedLinkFiltered.py /home/chrisgre/scripts/simulations/hg19.U87/simulation.${s}/oRNA.data.1.1.${c}.${i} /home/chrisgre/scripts/simulations/hg19.U87/simulation.${s}/all.aligned.ids.filtered.1.1.${c}.${i} /home/chrisgre/scripts/simulations/hg19.U87/simulation.${s}/oRNA.simSeqs
		done
	done
done
