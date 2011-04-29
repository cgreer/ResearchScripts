#!/bin/bash

for event in A3SS A5SS CE CI RI SE
#for event in A3SS A5SS CE CI
#for event in MXE MSE
do
    echo "Doing $event..."
    python2.4 classify_standard_event.py -a hg17 -d /r100/burge/shared/splice_graphs/hg17/FlatFiles -e $event -p SpliceGraphsCanonical_ -o StandardEventsCanonical_ -f /r100/burge/stadler/projects/splice_graphs/filter_loci.tsv
    python2.4 classify_standard_event.py -a hg17 -d /r100/burge/shared/splice_graphs/hg17/FlatFiles -e $event -p SpliceGraphsFiltered_ -o StandardEventsFiltered_ -f /r100/burge/stadler/projects/splice_graphs/filter_loci.tsv
    python2.4 classify_standard_event.py -a hg17 -d /r100/burge/shared/splice_graphs/hg17/FlatFiles -e $event -p SpliceGraphsTweaked_ -o StandardEventsTweaked_ -f /r100/burge/stadler/projects/splice_graphs/filter_loci.tsv

    #python2.4 classify_standard_event.py -a mm6 -d /r100/burge/shared/splice_graphs/mm6/FlatFiles -e $event -p SpliceGraphsCanonical_ -o StandardEventsCanonical_
    #python2.4 classify_standard_event.py -a mm6 -d /r100/burge/shared/splice_graphs/mm6/FlatFiles -e $event -p SpliceGraphsFiltered_ -o StandardEventsFiltered_
    #python2.4 classify_standard_event.py -a mm6 -d /r100/burge/shared/splice_graphs/mm6/FlatFiles -e $event -p SpliceGraphsTweaked_ -o StandardEventsTweaked_
done
