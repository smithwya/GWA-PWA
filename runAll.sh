#!/bin/bash
infile=2ch_kmat-nom_rhoN-nom
polefile=Poles/poletest
pWave=P
remin=0
remax=1
repts=100
immin=0
immax=1
impts=100
polecutoff=1
sheet=true
nJobs=1
numfits=15
runtime=6:00:00

mkdir -p Fits Plots Poles
make clean && make

for((i=1; i<=nJobs; i++))
do
sbatch --time=$runtime submit.script $infile $polefile $pWave $remin $remax $repts $immin $immax $impts $polecutoff $sheet "$i" $numfits
done
