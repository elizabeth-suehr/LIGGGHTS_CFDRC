#!/bin/bash
for i in 0 1 2 3 
do 
    qsub PBS_shear_rod6_${i}.script
   
done
