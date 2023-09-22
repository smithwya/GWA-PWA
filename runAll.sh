#!/bin/bash
nJobs=400
nFits=500
formatpath=Data/fitformat_exclusive.txt
fitpath=Fits/

mkdir -p Fits

make

for((i=1; i<=nJobs; i++))
do
sbatch submit.script "$i" $nFits $formatpath $fitpath
done
