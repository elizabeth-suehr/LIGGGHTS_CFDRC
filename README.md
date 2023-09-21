# LIGGGHTS_CFDRC
>  This code is a work in progress. Assume it only runs what has been tested and validated.

Results for the current dev branch should always be in [./scripts](./scripts/readme.md).
As of right now a push into the dev branch should indicated that Single Spheres, Rod2, Rod4, Rod6, and multisphere13 have been tested and verified. This list will eventually grow to include mixtures of single and multispheres, superquadric spheres, and superqudric cylinders.


# Current Progress
### Particle Types
|            | Single-Spheres | Multi-Spheres | Mixtures | SuperQuadrics |
|   -        | -------------  | ------------- | -------  | ------------- |
| **Starts** |     &#9745;    |     &#9745;   | &#9745;  |     &#9745;   |
| **Convergence** |     &#9745;    |     &#9745;   | &#9745;  |     &#9745;   |
| **Verified** |     &#9745;    |     &#9745; [complete](#multi_sphere)   | &#9744; [incomplete](#mixture)  |    &#9744;    |

### MPI Setting
> **Warning:**
> **Different MPI settings run very different code**. 

i.e. the current code should not run and will give printouts if you are using the LEBC with only 1 processor core in the y direction.  This isn't a very complicated fix, but if you don't have at least 2 processors cores to run this code, you probably don't need it.

The runs are now being done on a variety of MPI settings, and they all seem to be working. (execpt y=1)

<a name="multi_sphere"></a>
# Multi-Sphere
Rod2,4,6 and multisphere13 work is [here](./scripts/readme.md). For multispheres to be completed the following tests should compare well to literature.

- rod2  &#9745;
- rod4  &#9745;
- rod6  &#9745;
- multisphere13 (should give close to single-sphere results)  &#9745;


<a name="mixture"></a>
# Mixture
Curretly we have no validated results for mixtures because lebc.py isn't setup to generate mixtures input scripts correctly. 
The LIGGGHTS code should not need edited to get these results.

<a name="superquadrics"></a>
# Superquadrics
Curretly we have no verified results for superquadrics. Based on the few tests done with these shapes the collisonal stress tally is not working along with the fix_init's calcualtion of the overlap value.  The simulations seem to run and coverge to steady state values at least for translation stress.

## Notes
The LEBC Code was taken from fork of LIGGGHTS on github, the code on github has lots of bugs fixed by CFDRC and will not run most simulations we have tested in this repo.
https://github.com/SueHeir/LIGGGHTS-LEBC 

The FixInit and FixCool code was taken from a previous gitlab repo within CFDRC
[here](https://code.cfdrc.com/jason.howison/liggghts-custom).

