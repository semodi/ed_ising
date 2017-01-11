#!/bin/bash
#Before running this script, make sure that all the required files have been compiled

cp ed_config_readable.dat ed_config.dat

# Loop through symmetry representations (e.g. momenta); adjust upper limit
for (( i=1; i <= 1; i++ ))
do
 ./ed_update
 ./build_dict
 ./ising
done

