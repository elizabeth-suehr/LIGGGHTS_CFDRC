#!/bin/bash
for j in 0 1 2 3 4
do
    qsub PBS_shear_curl_3_${j}.script

done
