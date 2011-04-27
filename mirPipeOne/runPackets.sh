#!/bin/bash

for i in {1..10}
do
  python stageTWO.py -n $i &
  sleep 10
done
