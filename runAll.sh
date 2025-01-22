#!/bin/bash
infile=2ch_kmat-nom_rhoN-nom
nJobs=1
numfits=15
fitsfolder=/N/slate/giofoti/GWA-PWA/Fits/$(date +%F)/$infile
runtime=6:00:00
memory=1G
ncpus=4
mkdir -p $fitsfolder
make clean && make

for((i=1; i<=nJobs; i++))
do
sbatch --time=$runtime --mem=$memory --cpus-per-task=$ncpus -e /N/slate/giofoti/GWA-PWA/Fits/$(date +%F)/$infile-err.$i -o /N/slate/giofoti/GWA-PWA/Fits/$(date +%F)/$infile-out.$i submit.script "$infile" $i $numfits $fitsfolder $ncpus
done
