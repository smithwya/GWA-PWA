#!/bin/bash
# $1 job number
# $2 number of fits to perform
# $3 fit format file path
# $4 path to save valid fits

#runs jobs
for ((i=1; i<= $2; i++))
do
./bin/GWA $1 "$i" $3 $4
done
