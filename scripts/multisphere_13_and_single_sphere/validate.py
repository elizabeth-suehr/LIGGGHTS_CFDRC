#####################################################################################################
# Begin script, load the lebc function library
import sys
import os
sys.path.append("../.")
import lebc
import numpy as np
# #####################################################################################################
# # Auto Generate some inputs

# #Name of test is the directory file name, passed from validation.py

#Name of test is the directory file name, passed from validation.py
files = ["particle_data.txt",""]

particles_num_spheres = [13,1]



#####################################################################################################
# Test Specific inputs if needed


#####################################################################################################
# User inputs for LIGGGHTS simulations
nspecies        = 2					# Number of species to use
proc	 	= np.array([4,4,2])		 	# processor layout
ptype	 	= np.array(files)        	# Particle type, blank is sphere, file names means multisphere
# ptype	 	= np.array(["","","particle"])      	# Particle type, blank is sphere, file names means multisphere
nsph            = np.array(particles_num_spheres)                         # number of spheres for body, 1 for pure sphere (not body)
diam	 	= np.array([0,0.0001319*2])      		# Particle diameters, can be anything here for multisphere as it is computed before use
dens	 	= np.array([2500,2500])			# Particle densities
cres	 	= np.array([0.95, 0.95])		 	# Coefficients of restitution
amin            = 0.05					# Minimum volume fraction per species
amax         	= 0.5    	  			# Maximum volume fraction per species
asmax           = 0.68					# Maximum total volume fraction
distribution    = 'equal_ratio'     				# 'linear' or 'log' for alpha distribution 'equal_ratio' for each species takes up the same volume fraction as others species
npoints   	= 6					# Number of volume fractions, must be > 1
#nparticles      = 5000					# Total number of particles to use
filestr         = "multi_and_single" 	  			# Base filename string for input files
E		= 8.7e9					# LIGGGHTS parameter Young's modulus for all particles
nu		= 0.3					# LIGGGHTS parameter Poisson's ratio for all particles
steps     	= 15000000      			# Number of simulation steps after initialization
shearStr        = 100.0					# Shear strain rate to apply for shear simulations
lens            = 3.0e-4                               # representative length scale, bounding box length/largest size. Can be used to control #particles for a fixed volume
lebc_band       = 0.15
lebc_latt       = 0
gtemp           = 1.0                                   # granular temperature, used for cooling
rinit_xmax      = 1e-4                                  # initializer, maximum limit distance per timestep
rinit_reset     = 1000                                   # initializer, number of steps before resetting velocities to zero
rinit_th        = 1e-6                                  # initializer, threshold at which to exit the pre-run step (\sum{0.5*m_i*V_i^2}/\sum{m_i})
rinit_ns        = 1e5                                   # initializer, total maximum number of time steps to run the pre-run step
dtrelax         = 1.0                                   # time step relaxation (dt=dtrelax*dt_stable)
nbins_x         = 30                                    # number of particle col bins / xdir (shear is half, so it must be an even number). this is best if it is a multiple of the number of cores. 
# load inputs

lebc.set_input(nsph,proc,ptype,nspecies,diam,dens,cres,amin,amax,asmax,distribution,npoints,lens,lebc_band,lebc_latt,filestr,E,nu,steps,shearStr,gtemp,dtrelax,nbins_x,"true")
# Check for errors in input parameters
lebc.check_input()
# Calculations to set parameters specified in inputs
lebc.calculate_params()
# custom lebuf size (default is 0.001), Ly*lebuf is distance a particle must travel outside of the lebc before it is moved
lebc.set_lebuf(1e-3)
# # custom control of the initializer method, used in both cooling and shear sims
lebc.set_init_params(rinit_xmax,rinit_reset,rinit_th,rinit_ns)
# #************************************************************************************************************#

# BLOCK 1: Example to generate input files and qsub scripts to run all jobs
##################################################################
# Generate LIGGGHTS shear simulation input files
#lebc.shear_sim_files()
# # Generate LIGGGHTS cooling simulation input files
# # If a different processor layout is desired for cooling:
# proc = np.array([3,3,3])
# lebc.set_proc(proc)
# lebc.cooling_sim_files()
# ##################################################################
# # Generate PBS scripts to run each simulation
#lebc.pbs_files(1,6,6,"","16:00:00","normal","liz.suehr","lmp_auto","shear", "module load /home/liz.suehr/modules/openmpi/mvapich23-x86_64\nmodule load /home/liz.suehr/modules/liggghts/3.0\n")
# # lebc.pbs_files(1,28,27,"bro_ele","48:00:00","long","e1530","lmp_autoP19","cooling")

# #************************************************************************************************************#

# # BLOCK2: Example to post-process simulation data into GGFS compatible file
# # Load the shear simulation data to memory
lebc.load_shear_data()
# # # Load the cooling simulation data to memory
# # lebc.load_cooling_data()
# # plot data to check (compare to mixture-averaged density/diameter kinetic theory for sphere)
lebc.plot_shear_data()
# # lebc.plot_cooling_data()
# # To check if cooling simulations are actually in the HCS regime
for i in range(0,npoints):
    lebc.plot_shear_case_data(i)
# # # to visually inspect data before writing the HDF5 file
# # lebc.inspect_data()
# # # # To write the ggfs hdf5 file, both cooling and shear data needs to be loaded
# # lebc.data_to_ggfs("testRun")

# #************************************************************************************************************#
