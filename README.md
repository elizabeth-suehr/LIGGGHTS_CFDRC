# LIGGGHTS_CFDRC
LIGGGHTS documentation can be found at https://www.cfdem.com/media/DEM/docu/Section_intro.html

The new Lees-Edwards boundary condition command is
```
fix             leboundary all lebc 100.0 true gtemp 1e-09 ave_reset 60000 body_data body_data_name 10000
```
where lebc 100.0 is the shear strain rate, gtemp is contant to control the inital randomize velocity ( helps with low volume fraction runs reaching steady state faster), ave_reset is how many time steps the stress is averaged over, body_data is how often you print out all the body data


Initalization at high volume fractions is typically difficult, in this code base the insert/pack command as the option 'overlapcheck yes' disabled. This is because we have a relaxation method to remove overlaps with the following command
```
fix             limcheck all slowinit xmax 1.00e-07 reset 1000 threshold 1.00e-07 start_dt 1.000e-15 end_dt 5.542e-09
```

This is ran and then removed before the LEBC fix is added. For many examples of using this codebase please see https://github.com/elizabeth-suehr/dem-scripts


# Current Progress
### Particle Types
|            | Single-Spheres | Multi-Spheres | Mixtures | SuperQuadrics |
|   -        | -------------  | ------------- | -------  | ------------- |
| **Starts** |     &#9745;    |     &#9745;   | &#9745;  |     &#9745;   |
| **Convergence** |     &#9745;    |     &#9745;   | &#9745;  |     &#9745;   |
| **Validated** |     &#9745;    |      &#9745;   | [incomplete](#mixture)  |    &#9744;    |





