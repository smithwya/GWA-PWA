#!/bin/bash
nJobs=100
nFits=100
formatpath=Data/fitformat_inclusive.txt
fitpath=Fits/

mkdir -p Fits

make

for((i=1; i<=nJobs; i++))
do
sbatch submit.script "$i" $nFits $formatpath $fitpath
done
