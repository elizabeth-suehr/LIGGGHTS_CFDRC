# LIGGGHTS_CFDRC
>  This code is a work in progress. Assume it only runs what has been tested and validated.

Results for the current dev branch should always be in [./scripts](./scripts/readme.md).
As of right now a push into the dev branch should indicated that Single Spheres, Rod2, and multisphere13 have been tested and validated. This list will eventually grow to include mixtures of single and multispheres, superquadric spheres.


# Current Progress
### Particle Types
|            | Single-Spheres | Multi-Spheres | Mixtures | SuperQuadrics |
|   -        | -------------  | ------------- | -------  | ------------- |
| **Starts** |     &#9745;    |     &#9745;   | &#9745;  |     &#9745;   |
| **Convergence** |     &#9745;    |     &#9745;   | &#9745;  |     &#9745;   |
| **Validated** |     &#9745;    |     [incomplete](#multi_sphere)   | [incomplete](#mixture)  |    &#9744;    |

### MPI Setting
> **Warning:**
> **Different MPI settings run very different code**. 

i.e. the current code should not run and will give printouts if you are using the LEBC with only 1 processor core in the y direction.  This isn't a very complicated fix, but if you don't have at least 2 processors cores to run this code, you probably don't need it.

Currently the single-sphere and multi-sphere runs are done on a MPI settings of 3x2x1 (x y z).  It is likely that any MPI x setting, z setting, and y setting greater than 1 should all work. That being said, if it hasn't been validated then that specific MPI setting might not work. 



<a name="multi_sphere"></a>
# Multi-Sphere
Rod2 and multisphere13 work is [here](./scripts/readme.md). For multispheres to be completed the following tests should compare well to literature.

- rod2  &#9745;
- rod4  &#9744;
- rod6  &#9744;
- multisphere13 (should give close to single-sphere results)  &#9745;


<a name="mixture"></a>
# Mixture
Curretly we have no validated results for mixtures because lebc.py isn't setup to generate mixtures input scripts correctly. 
The LIGGGHTS code should not need edited to get these results.


## Notes
The LEBC Code was taken from fork of LIGGGHTS on github, the code on github has lots of bugs fixed by CFDRC and will not run most simulations we have tested in this repo.
https://github.com/SueHeir/LIGGGHTS-LEBC 

The FixInit and FixCool code was taken from a previous gitlab repo within CFDRC
[here](https://code.cfdrc.com/jason.howison/liggghts-custom).

