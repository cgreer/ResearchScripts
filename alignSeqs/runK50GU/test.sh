#!/bin/bash

gawk '$16=="T"' $1 > $2
