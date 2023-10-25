#!/bin/bash
nJobs=400
nFits=25
formatpath=Data/fitformat_inclusive.txt
fitpath=Fits/
runtime=24:00:00

mkdir -p Fits Plots

make

for((i=1; i<=nJobs; i++))
do
sbatch --time=$runtime submit.script "$i" $nFits $formatpath $fitpath
done
