#!/bin/bash

mkdir output/run1
mkdir output/run1/u
for i in {0..2}  
do
   ./phosim input/run1_3650/lightcurve_${i}.icat -c input/quickbackground -o output/run1/u -e 0 -i lsst
   echo "Generated image for time step $i"
done

