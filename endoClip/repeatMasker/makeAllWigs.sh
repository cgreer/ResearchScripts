#!/bin/bash

for i in {1..22}; do $HOME/exec/qJobX9.sh $HOME/exec/qDo.sh makeRepeatWigs.py makeRepeatWigs chr$i 1 repeatMasker.txt repeatWigs; done 
for i in {1..22}; do $HOME/exec/qJobX9.sh $HOME/exec/qDo.sh makeRepeatWigs.py makeRepeatWigs chr$i -1 repeatMasker.txt repeatWigs; done 
for i in X Y M; do $HOME/exec/qJobX9.sh $HOME/exec/qDo.sh makeRepeatWigs.py makeRepeatWigs chr$i 1 repeatMasker.txt repeatWigs; done 
for i in X Y M; do $HOME/exec/qJobX9.sh $HOME/exec/qDo.sh makeRepeatWigs.py makeRepeatWigs chr$i -1 repeatMasker.txt repeatWigs; done 
