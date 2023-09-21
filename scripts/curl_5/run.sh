#!/bin/bash
for j in 0 1 2 3 4 5 6
do
    qsub PBS_shear_curl_5_${j}.script

done
