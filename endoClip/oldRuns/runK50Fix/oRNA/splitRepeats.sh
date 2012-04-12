#!/bin/bash


for j in 50 60 70 80; do for i in {0..5}; do gawk '$22=="T" || $23<36' oRNA.data.1.1.$j.$i > oRNA.data.1.1.$j.$i.repeat; done; done
for j in 50 60 70 80; do for i in {0..5}; do gawk '$22=="F" && $23>35' oRNA.data.1.1.$j.$i > oRNA.data.1.1.$j.$i.norepeat; done; done
