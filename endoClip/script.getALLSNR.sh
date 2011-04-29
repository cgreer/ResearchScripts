#!/bin/bash

for i in 50 60; do for k in {0..3}; do for j in norepeat repeat; do python script.getSNRStats.py runK50GU/oRNA.data.1.1.$i.$k.$j >> ALL.SNR; done; done; done
