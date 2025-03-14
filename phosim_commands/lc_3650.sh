#!/bin/bash

for i in {0..124}  
do
   ./phosim custom/lsst/3_14_3650/lightcurve_${i}.icat -c custom/lsst/quickbackground -o custom/lsst/3_14_3650_output -e 0 -i lsst
   echo "Generated image for time step $i"
done

