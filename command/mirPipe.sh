#!/bin/bash

echo 'Defining Clusters'
cd /home/chrisgre/scripts/defineClusters
python defineClusters.py $1

echo 'Filtering Out Known MicroRNAs/snoRNAs'
cd /home/chrisgre/scripts/filterKnown
python filterOutKnown.py $1

echo 'Updating Results File/Running Initial Sort'
cd /home/chrisgre/scripts/updateFPF
python addPeriods.py $1
python updateDensity.py $1
python updateOverlaps.py $1
python updateClusterInfo.py $1
python sortResults.py $1

echo 'Splitting Results into Exons and Introns'
cd /home/chrisgre/scripts/splitExonIntron
python splitExonsIntrons.py $1

echo 'Finding Highest Read Count Per Cluster'
cd /home/chrisgre/scripts/readDensity

echo '  Exons...'
python getMaxReadDensity.py E  $1
echo '  Introns...'
python getMaxReadDensity.py I $1


echo 'Collecting Noise Data'
cd /home/chrisgre/scripts/noisyRanges
echo '  Exons...'
python exonNoisy.py $1
echo '  Introns...'
python intronNoisy.py $1

echo 'Updating PVals For Prediction File'
python updateSignalNoise.py E $1
python updateSignalNoise.py I $1

echo 'Sorting (FINAL)'
cd /home/chrisgre/scripts/updateFPF
python finalSort.py E $1
python finalSort.py I $1

echo 'Locating Predictions w/ Canonical Peaks'
cd /home/chrisgre/scripts/seqPeaks
python bestSinglePeakPlus.py E $1
python bestSinglePeakPlus.py I $1

echo 'Done'
