cd ../src
# make clean-all
make auto -j 6
cd ../test
rm lmp_auto
rm restart_forced_*
ln -s ../src/lmp_auto lmp_auto
#./lmp_auto < in.single_sphere_1_1_1 
#mpiexec -np 2 lmp_auto < in.single_sphere_1_2_1
#mpiexec -np 2 lmp_auto < in.single_sphere_2_1_1
#mpiexec -np 3 lmp_auto < in.single_sphere_3_1_1

#mpiexec -np 6 lmp_auto < in.single_sphere_3_2_1
#mpiexec -np 4 lmp_auto < in.single_sphere_4_1_1


#./lmp_auto < in.multi_sphere_1_1_1 
#mpiexec -np 2 lmp_auto < in.multi_sphere_1_1_1



#mpirun -n 6 xterm -e gdb -x commands.gdb ./lmp_auto
#mpiexec -np 4 lmp_auto < in.multi_sphere_4_1_1

#mpiexec -np 6 lmp_auto < in.single_sphere_2_2_1 > stdout_ss_2_2_1.txt
mpiexec -np 6 lmp_auto < in.multi_sphere_2_2_1 > stdout_ms_2_2_1.txt


