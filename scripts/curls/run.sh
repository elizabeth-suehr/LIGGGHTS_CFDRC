#!/bin/bash
for i in 0 2 3 4 
do 
    #!/bin/bash
    for j in 0 1 2 3 4 5
    do
        qsub PBS_shear_curl_${i}_${j}.script

    done
   
done
