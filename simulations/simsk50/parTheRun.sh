#!/bin/bash

simDir=$1

mkdir $simDir/pRuns
for i in {0..9}; do mkdir ${simDir}/pRuns/run.0$i; done
/home/chrisgre/exec/splitFile.sh ${simDir}/all.aligned.nodups.truncated ${simDir} 10
for i in {0..9}; do cp ${simDir}/all.aligned.nodups.truncated.splitted.0$i ${simDir}/pRuns/run.0$i/all.aligned.nodups.truncated; done
for i in ${simDir}/pRuns/run.0*; do mkdir $i/oRNA; done
for i in ${simDir}/pRuns/run.0*; do mkdir $i/aDir; done
for i in ${simDir}/all.aligned.nodups.splitted*; do rm $i; done
for i in {0..9}; do cp ${simDir}/oRNA.simSeqs ${simDir}/pRuns/run.0$i/oRNA.simSeqs; done

