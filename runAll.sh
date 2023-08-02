#!/bin/bash
nJobs=100
nFits=30
formatpath=Data/fitformat_exclusive.txt
fitpath=Fits/

make

for((i=1; i<=nJobs; i++))
do
sbatch submit.script "$i" $nFits $formatpath $fitpath
done
