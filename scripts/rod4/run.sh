#!/bin/bash
for i in 0 1 2 3 4 5 6 7 
do 
    qsub PBS_shear_rod4_${i}.script
   
done
