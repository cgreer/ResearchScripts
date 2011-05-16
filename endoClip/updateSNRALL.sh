#!/bin/bash

for i in 50 60; do for j in {0..5}; do python updateSNRFlat.py runK1/oRNA.data.1.1.$i.$j; done; done
