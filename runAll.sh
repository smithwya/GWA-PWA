#!/bin/bash
nJobs=30
nFits=10
formatpath=Data/fitformat.txt
fitpath=Fits/

make

for((i=1; i<=nJobs; i++))
do
sbatch submit.script "$i" $nFits $formatpath $fitpath
done
