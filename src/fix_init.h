/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    This file is from LAMMPS
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(slowinit,FixELimit)

#else

#ifndef LMP_FIX_INIT_H
#define LMP_FIX_INIT_H

#include "fix.h"
#include <vector>

namespace LAMMPS_NS {

class FixELimit : public Fix {
 public:
  FixELimit(class LAMMPS *, int, char **);
  ~FixELimit();
  int setmask();
  void init();
  void initial_integrate(int);
  void final_integrate();
  void end_of_step();
  void reset_dt();
  void reset_vel();
  double check_ke();
  double compute_scalar();
  void post_force(int dummy);

  double variableCOR(double overlap);
 private:
  double global_minOverlap;
  double start_dt = 1e-16;
  double end_dt = 1e-10;

  double max_cor = 0.1;
  double cke_mass[10];
  int cke_icheck;
  std::vector<double> dtii;
  std::vector<double> cor0;
  double dtv,dtf;
  double *step_respa;
  int mass_require,ncount;
  double xlimit,vlimitsq;

  int check_int,check_n_int;
  int relflag; 
  double threshold;
  class FixMultisphere *fix_ms;
};

}

#endif
#endif
