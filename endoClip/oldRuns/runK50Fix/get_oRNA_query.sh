#!/bin/bash

gawk '$9>1.2 && $11<7 && $12<7 && $13=="F" {print $1, $4}' $1 | tr ' ' '\t'
