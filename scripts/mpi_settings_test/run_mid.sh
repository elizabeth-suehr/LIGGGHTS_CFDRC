#!/bin/bash


for j in 221
do
    for i in 13 43 93 163
    do
        qsub PBS_shear_multisphere${i}_MPI${j}_0.script
    done
    qsub PBS_shear_s_sphere_MPI${j}_0.script

done


for j in 221
do
    for i in 13 43 93 163
    do
        qsub PBS_shear_multisphere${i}_MPI${j}_1.script
    done
    qsub PBS_shear_s_sphere_MPI${j}_1.script

done


# for j in 884
# do
#     for i in 13 43 93 163
#     do
#         qsub PBS_shear_multisphere${i}_MPI${j}_0.script
#     done
#     qsub PBS_shear_s_sphere_MPI${j}_0.script

# done


# for j in 884
# do
#     for i in 13 43 93 163
#     do
#         qsub PBS_shear_multisphere${i}_MPI${j}_1.script
#     done
#     qsub PBS_shear_s_sphere_MPI${j}_1.script

# done
