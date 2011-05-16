#!/bin/bash

simDir=$1

mkdir $simDir/pRuns2
for i in {0..49}; do mkdir ${simDir}/pRuns2/run.$i; done
/home/chrisgre/exec/splitFile.sh ${simDir}/all.aligned.0.5.apr14.ids ${simDir} 50
#rename the beginning ten
for i in {0..9}; do mv ${simDir}/all.aligned.0.5.apr14.ids.splitted.0$i ${simDir}/all.aligned.0.5.apr14.ids.splitted.$i; done
for i in {0..49}; do mv ${simDir}/all.aligned.0.5.apr14.ids.splitted.$i ${simDir}/pRuns2/run.$i/all.aligned.nodups.truncated; done
for i in ${simDir}/pRuns2/run.*; do mkdir $i/oRNA; done
for i in ${simDir}/pRuns2/run.*; do mkdir $i/aDir; done
for i in {0..49}; do cp ${simDir}/oRNA.simSeqs ${simDir}/pRuns2/run.$i/oRNA.simSeqs; done

