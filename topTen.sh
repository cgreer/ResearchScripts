#!/bin/bash

echo

for i in compute-7-1 compute-7-2 compute-7-3 compute-7-4
do

echo $i
ssh $i ps -eo pcpu,pmem,pid,user,args | sort -k 1 -r | head -12
echo

done
