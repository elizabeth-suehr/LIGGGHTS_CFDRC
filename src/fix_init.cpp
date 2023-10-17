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
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fix_init.h"
#include "atom.h"
#include "force.h"
#include "pair.h"
#include "update.h"
#include "respa.h"
#include "modify.h"
#include "comm.h"
#include "error.h"
#include "fix_multisphere.h"
#include "fix_template_multisphere.h"
#include "particleToInsert_multisphere.h"
#include "multisphere_parallel.h"
#include "signal_handling.h"
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace LAMMPS_NS;
using namespace FixConst;
static const char *COEFFICIENT_RESTITUTION = "coefficientRestitution";
/* ---------------------------------------------------------------------- */
#define maxNsteps 20000
#define MIN_OVERLAP 0.6

FixELimit::FixELimit(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  // if (comm->me==0)std::cout << "pre-run, create \n";
  // if (narg != 5) error->all(FLERR,"Illegal fix nve/limit command");

  time_integrate = 1;
  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  check_int = 10;
  int iarg = 3;
  xlimit = 0.1;
  if (strcmp(arg[iarg], "xmax") == 0)
  {
    if (narg < iarg + 1)
      error->fix_error(FLERR, this, "Not enough arguments for 'xmax' in Limit \n");
    xlimit = force->numeric(FLERR, arg[iarg + 1]);
    iarg += 2;
    relflag = 0;
  }

  if (strcmp(arg[iarg], "reset") == 0)
  {
    if (narg < iarg + 1)
      error->fix_error(FLERR, this, "Not enough arguments for 'reset' in Limit \n");
    check_int = force->numeric(FLERR, arg[iarg + 1]);
    iarg += 2;
  }

  threshold = 1e-6;
  if (strcmp(arg[iarg], "threshold") == 0)
  {
    if (narg < iarg + 1)
      error->fix_error(FLERR, this, "Not enough arguments for 'threshold' in Limit \n");
    threshold = force->numeric(FLERR, arg[iarg + 1]);
    iarg += 2;
  }

  if (strcmp(arg[iarg], "start_dt") == 0)
  {
    if (narg < iarg + 1)
      error->fix_error(FLERR, this, "Not enough arguments for 'start_dt' in Limit \n");
    start_dt = force->numeric(FLERR, arg[iarg + 1]);
    update->dt = start_dt; // * 1e5;
    update->timestep_set = true;
    iarg += 2;
    relflag = 0;
  }

  if (strcmp(arg[iarg], "end_dt") == 0)
  {
    if (narg < iarg + 1)
      error->fix_error(FLERR, this, "Not enough arguments for 'end_dt' in Limit \n");
    end_dt = force->numeric(FLERR, arg[iarg + 1]);
    iarg += 2;
    relflag = 0;
  }

  ncount = 0;

  fix_ms = NULL;
  // if (comm->me==0)std::cout << "pre-run, create done\n";
}

FixELimit::~FixELimit()
{
  if (fix_ms)
  {
    fix_ms->force_volume_fraction = false;

    std::vector<double> xcm;
    std::vector<double> quat;

    double temp[3] = {0.0, 0.0, 0.0};
    double temp_quat[4] = {0.0, 0.0, 0.0, 0.0};
    for (int i = 0; i < fix_ms->n_body(); i++)
    {
      fix_ms->data().xcm(temp, i);
      fix_ms->data().quat(temp_quat, i);
      xcm.push_back(temp[0]);
      xcm.push_back(temp[1]);
      xcm.push_back(temp[2]);
      quat.push_back(temp_quat[0]);
      quat.push_back(temp_quat[1]);
      quat.push_back(temp_quat[2]);
      quat.push_back(temp_quat[3]);
    }

    std::cout << xcm[0] << quat[0];

    std::vector<int> counts_recv_xcm(comm->nprocs);
    std::vector<int> displacements_xcm(comm->nprocs);

    std::vector<int> counts_recv_quat(comm->nprocs);
    std::vector<int> displacements_quat(comm->nprocs);

    int size_xcm = xcm.size();
    int size_quat = quat.size();

    MPI_Gather(&size_xcm,
               1,
               MPI_INT,
               &counts_recv_xcm[0],
               1,
               MPI_INT,
               0,
               world);

    MPI_Gather(&size_quat,
               1,
               MPI_INT,
               &counts_recv_quat[0],
               1,
               MPI_INT,
               0,
               world);

    int displs = 0;
    for (int i = 0; i < comm->nprocs; ++i)
    {
      displacements_xcm[i] = displs;
      displs += counts_recv_xcm[i];
    }

    std::vector<double> xcm_all(displs);

    int displs_quat = 0;
    for (int i = 0; i < comm->nprocs; ++i)
    {
      displacements_quat[i] = displs_quat;
      displs_quat += counts_recv_quat[i];
    }
    std::vector<double> quat_all(displs_quat);

    MPI_Gatherv(&xcm[0],
                xcm.size(),
                MPI_DOUBLE,
                &xcm_all[0],
                &counts_recv_xcm[0],
                &displacements_xcm[0],
                MPI_DOUBLE,
                0,
                world);

    MPI_Gatherv(&quat[0],
                quat.size(),
                MPI_DOUBLE,
                &quat_all[0],
                &counts_recv_quat[0],
                &displacements_quat[0],
                MPI_DOUBLE,
                0,
                world);

    if (comm->me == 0)
    {
      std::ofstream ofs;

      std::string filename = "xcm_quat";
      filename += std::to_string(fix_ms->n_body_all());
      filename += ".txt";
      ofs.open(filename.c_str(), std::ofstream::out | std::ofstream::trunc);

      for (int i = 0; i < fix_ms->n_body_all(); i++)
      {
        // std::cout << i + 1 << " " << xcm_all[i * 3] << " " << xcm_all[i * 3 + 1] << " " << xcm_all[i * 3 + 2] << " " << quat_all[i * 4] << " " << quat_all[i * 4 + 1] << " " << quat_all[i * 4 + 2] << " " << quat_all[i * 4 + 3] << "\n";
        ofs << i + 1 << " " << xcm_all[i * 3] << " " << xcm_all[i * 3 + 1] << " " << xcm_all[i * 3 + 2] << " " << quat_all[i * 4] << " " << quat_all[i * 4 + 1] << " " << quat_all[i * 4 + 2] << " " << quat_all[i * 4 + 3] << "\n";
      }
      ofs.close();
    }
  }

  SignalHandler::request_exit = false;
}

/* ---------------------------------------------------------------------- */

int FixELimit::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= END_OF_STEP;
  return mask;
}

void FixELimit::post_force(int dummy)
{
  int *mask = atom->mask;
  int *itype = atom->type;
  double overlap = {1.0};
  if (force->pair)
  {
    double **fl_overlap = force->pair->overlap;
    int npair = atom->nlocal;
    for (int i = 0; i < npair; i++)
      overlap = std::min(overlap, fl_overlap[i][0]);
  }
  // if (update->ntimestep % 1000 == 0)
  //   for (int i=0; i<comm->nprocs; i++){
  //     if (i==comm->me)
  //       std::cout<<i<<", overL="<<overlap<<"\n";
  //     MPI_Barrier(MPI_COMM_WORLD);
  //   }

  MPI_Allreduce(&overlap, &global_minOverlap, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
}
/* ---------------------------------------------------------------------- */

void logSpace(double a, double b, int np, double *ISO)
{
  double minISO = log(a);
  double maxISO = log(b);
  double dx = (maxISO - minISO) / (np - 1.);

  ISO[np - 1] = maxISO;
  for (int i = 0; i < np - 1; i++)
    ISO[np - (i + 2)] = ISO[np - (i + 1)] - dx;
  for (int i = 0; i < np; i++)
    ISO[i] = exp(ISO[i]);
}

void FixELimit::init()
{
  global_minOverlap = 0.5; // just to get it started
  cke_icheck = 2;
  check_n_int = 0;
  reset_vel();
  int nnstep = maxNsteps;
  std::vector<double> dt_slow(nnstep);
  logSpace(start_dt, end_dt * 9.9, nnstep, &dt_slow[0]);
  update->dt = start_dt;
  update->timestep_set = true;
  for (int i = 0; i < nnstep; i++)
    dtii.push_back(dt_slow[i]);
  // if (!fix_ms)
  fix_ms = static_cast<FixMultisphere *>(modify->find_fix_style("multisphere", 0));

  if (fix_ms)
  {
    fix_ms->force_volume_fraction = true;
  }

  // std::cout << "pre-run, init() \n";
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  vlimitsq = (xlimit / dtv) * (xlimit / dtv);
  ncount = 0;

  if (relflag == 1 && (!atom->radius_flag || !atom->rmass_flag))
    error->fix_error(FLERR, this, "using 'radius_ratio' needs per-atom radius and mass");

  if (strstr(update->integrate_style, "respa"))
    step_respa = ((Respa *)update->integrate)->step;
  // warn if using fix shake, which will lead to invalid constraint forces

  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style, "shake") == 0)
    {
      if (comm->me == 0)
        error->warning(FLERR, "Should not use fix nve/limit with fix shake");
    }
  force->pair->lebc_setup(); // MAYBE?
  force->pair->pre_init_lebc();
  if (fix_ms)
    fix_ms->lebc_highest_body_id = fix_ms->tag_max_body();
  // if (fix_ms) {
  //   // force->pair->lebc_setup_multisphere();
  //   // fix_ms->lebc_highest_body_id = fix_ms->tag_max_body();
  // }
  // // compute xlimit based on given time-step and size of particles
  // double dtm = update->dt;
  // double pdiam_min = 1.0e10;
  // double *radius = atom->radius;
  // int *mask = atom->mask;
  // int nlocal = atom->nlocal;
  // for (int i = 0; i < nlocal; i++)
  //   if (mask[i] & groupbit)
  //     pdiam_min = std::min(pdiam_min,radius[i]);
  // xlimit = pdiam_min/dtm/100.0;
  // if (comm->me==0) std::cout << "pre-run, init(): vlimit= "<< pdiam_min/dtm/100.0 <<" \n";
  {
    std::string tmp = "FixELimit";
    const int max_type = force->registry.max_type();
    MatrixProperty *coeffRestProp = force->registry.getMatrixProperty(COEFFICIENT_RESTITUTION, tmp.c_str());
    double **coeffRest = coeffRestProp->data;
    int id = 0;
    for (int i = 1; i < max_type + 1; i++)
    {
      for (int j = 1; j < max_type + 1; j++)
      {
        cor0.push_back(coeffRest[i][j]);
        coeffRest[i][j] = 0.1;
      }
    }
  }
}

/////////////////
double FixELimit::variableCOR(double overlap)
{

  max_cor = std::max(0.9 / (1.0 + exp(-17.5 * (overlap - 0.8))) + 0.1, max_cor);

  return max_cor;

  if (overlap < MIN_OVERLAP)
    return 0.1;
  else
  {
    double temp = std::min((0.95 - 0.1) / (0.95 - MIN_OVERLAP) * (overlap - MIN_OVERLAP) + 0.1, 0.99);

    return temp * (overlap * overlap * overlap);
  }
}
/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixELimit::initial_integrate(int vflag)
{
  // variable COR, from 0.1 to 0.95 to make sure minimal overlap is reached fast.
  {
    std::string tmp = "FixELimit";
    const int max_type = force->registry.max_type();
    MatrixProperty *coeffRestProp = force->registry.getMatrixProperty(COEFFICIENT_RESTITUTION, tmp.c_str());
    double **coeffRest = coeffRestProp->data;
    int id = 0;
    for (int i = 1; i < max_type + 1; i++)
    {
      for (int j = 1; j < max_type + 1; j++)
      {
        coeffRest[i][j] = variableCOR(global_minOverlap);
      }
    }
  }
  // if (comm->me == 0) std::cout << "pre-run, initial_integrate() \n";
  double dtfm, vsq, scale, rsq;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  if (rmass)
  {
    if (relflag == 1)
    {
      for (int i = 0; i < nlocal; i++)
      {
        if (mask[i] & groupbit)
        {
          dtfm = dtf / rmass[i];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];

          vsq = v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];
          rsq = radius[i] * radius[i];
          if (vsq > vlimitsq * rsq)
          {
            ncount++;
            scale = sqrt(vlimitsq * rsq / vsq);
            v[i][0] *= scale;
            v[i][1] *= scale;
            v[i][2] *= scale;
          }

          x[i][0] += dtv * v[i][0];
          x[i][1] += dtv * v[i][1];
          x[i][2] += dtv * v[i][2];
        }
      }
    }
    else
    {
      for (int i = 0; i < nlocal; i++)
      {
        if (mask[i] & groupbit)
        {
          dtfm = dtf / rmass[i];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];

          vsq = v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];
          if (vsq > vlimitsq)
          {
            ncount++;
            scale = sqrt(vlimitsq / vsq);
            v[i][0] *= scale;
            v[i][1] *= scale;
            v[i][2] *= scale;
          }

          x[i][0] += dtv * v[i][0];
          x[i][1] += dtv * v[i][1];
          x[i][2] += dtv * v[i][2];
        }
      }
    }
  }
  else
  {
    for (int i = 0; i < nlocal; i++)
    {
      if (mask[i] & groupbit)
      {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];

        vsq = v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];
        if (vsq > vlimitsq)
        {
          ncount++;
          scale = sqrt(vlimitsq / vsq);
          v[i][0] *= scale;
          v[i][1] *= scale;
          v[i][2] *= scale;
        }

        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixELimit::final_integrate()
{
  // if (comm->me == 0) std::cout << "pre-run, final_integrate() \n";
  double dtfm, vsq, scale, rsq;

  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  if (rmass)
  {
    if (relflag == 1)
    {
      for (int i = 0; i < nlocal; i++)
      {
        if (mask[i] & groupbit)
        {
          dtfm = dtf / rmass[i];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];

          vsq = v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];
          rsq = radius[i] * radius[i];
          if (vsq > vlimitsq * rsq)
          {
            ncount++;
            scale = sqrt(vlimitsq * rsq / vsq);
            v[i][0] *= scale;
            v[i][1] *= scale;
            v[i][2] *= scale;
          }
        }
      }
    }
    else
    {
      for (int i = 0; i < nlocal; i++)
      {
        if (mask[i] & groupbit)
        {
          dtfm = dtf / rmass[i];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];

          vsq = v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];
          if (vsq > vlimitsq)
          {
            ncount++;
            scale = sqrt(vlimitsq / vsq);
            v[i][0] *= scale;
            v[i][1] *= scale;
            v[i][2] *= scale;
          }
        }
      }
    }
  }
  else
  {
    for (int i = 0; i < nlocal; i++)
    {
      if (mask[i] & groupbit)
      {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];

        vsq = v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];
        if (vsq > vlimitsq)
        {
          ncount++;
          scale = sqrt(vlimitsq / vsq);
          v[i][0] *= scale;
          v[i][1] *= scale;
          v[i][2] *= scale;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixELimit::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  vlimitsq = (xlimit / dtv) * (xlimit / dtv);
}

/* ----------------------------------------------------------------------
   energy of indenter interaction
------------------------------------------------------------------------- */

double FixELimit::compute_scalar()
{
  double one = ncount;
  double all;
  MPI_Allreduce(&one, &all, 1, MPI_DOUBLE, MPI_SUM, world);
  return all;
}

/* ----------------------------------------------------------------------
   check kinetic energy
------------------------------------------------------------------------- */

double FixELimit::check_ke()
{
  // return 0;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double ke = 0.0;

  if (rmass)
  {
    for (int i = 0; i < nlocal; i++)
    {
      if (mask[i] & groupbit && (!fix_ms || fix_ms->belongs_to(i) < 0))
      {
        const double mult = 0.0; // halfstep ? update->dt*0.5/rmass[i] : 0.0;
        const double v0 = v[i][0] + mult * f[i][0];
        const double v1 = v[i][1] + mult * f[i][1];
        const double v2 = v[i][2] + mult * f[i][2];
        ke += rmass[i] * (v0 * v0 + v1 * v1 + v2 * v2);
      }
    }
  }
  else
  {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        ke += mass[type[i]] *
              (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);
  }
  double scalar = 0.0;
  MPI_Allreduce(&ke, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);
  scalar *= 0.5; // pfactor;

  // for multispheres we need to get the kinetic energy from the multisphere fix
  if (fix_ms)
    scalar += fix_ms->extract_ke();

  return scalar;
}

void FixELimit::reset_vel()
{
  double v_cm[3] = {0.0, 0.0, 0.0};
  int natom = atom->nlocal + atom->nghost;
  if (fix_ms)
  {
    for (int i = 0; i < fix_ms->n_body(); i++)
    {
      fix_ms->data().set_v_body(i, v_cm);
      double cm[3] = {0, 0, 0};
      fix_ms->data().reset_forces(0);
      fix_ms->data().set_angmom_via_omega_body(i, cm);
    }
    for (int i = 0; i < natom; i++)
    {
      if (fix_ms->belongs_to(i) == -1)
      {
        atom->v[i][0] = v_cm[0];
        atom->v[i][1] = v_cm[1];
        atom->v[i][2] = v_cm[2];
      }
    }
  }
  else
  {
    for (int i = 0; i < natom; i++)
    {
      atom->v[i][0] = v_cm[0];
      atom->v[i][1] = v_cm[1];
      atom->v[i][2] = v_cm[2];
    }
  }

  if (fix_ms)
  { // needed
    fix_ms->set_xv();
    fix_ms->out_rev_comm_x_v_omega();
    fix_ms->add_body_finalize();
  }
  // update->dt = update->dt*1e5;
  // update->timestep_set = true;
}

void FixELimit::end_of_step()
{

  if (update->ntimestep < maxNsteps)
  {
    double temp = dtii[update->ntimestep];
    update->dt = std::min(temp, end_dt * (10.0 * variableCOR(global_minOverlap)));
    if (fix_ms)
    {
      fix_ms->update_dt();
    }
  }
  else
  {
    update->dt = end_dt * (10.0 * variableCOR(global_minOverlap));

    if (fix_ms)
    {
      fix_ms->update_dt();
    }
  }

  reset_dt();
  // if (comm->me == 0) std::cout << "pre-run, step: "<< update->ntimestep<<", dt="<<update->dt<<", natoms="<<atom->natoms<<"\n";
  int check_exit = 0;
  // if (update->ntimestep % 1000 == 0)
  // {
  double cke = check_ke();
  // divide by mass to get a sense of relative velocities
  double mass_sum = 0.0;
  int *mask = atom->mask;
  int *itype = atom->type;
  int nlocal = atom->nlocal;

  if (fix_ms)
  {

    if (atom->rmass)
      for (int i = 0; i < atom->nlocal; i++)
      {
        if (mask[i] & groupbit && (atom->tag[i] <= atom->lebc_highest_atom_id && ((fix_ms) ? (fix_ms->belongs_to(i) == -1) : 1)))
        {
          double massp = atom->rmass[i];
          mass_sum += massp;
        }
      }
    else
      for (int i = 0; i < atom->nlocal; i++)
      {
        if (mask[i] & groupbit)
        {
          double massp = atom->mass[itype[i]];
          mass_sum += massp;
        }
      }

    for (int i = 0; i < fix_ms->n_body(); i++)
    {
      if (fix_ms->data().tag(i) != -1)
      {
        double massB = fix_ms->data().mass(i);
        mass_sum += massB;
      }
    }
  }
  else
  {
    if (atom->rmass)
      for (int i = 0; i < atom->nlocal; i++)
      {
        if (mask[i] & groupbit)
        {
          double massp = atom->rmass[i];
          mass_sum += massp;
        }
      }
    else
      for (int i = 0; i < atom->nlocal; i++)
      {
        if (mask[i] & groupbit)
        {
          double massp = atom->mass[itype[i]];
          mass_sum += massp;
        }
      }
  }

  double global_sum_mass;
  MPI_Allreduce(&mass_sum, &global_sum_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (update->ntimestep % 1000 == 0)
  {

    int id = 0;
    double coravg = 0.0;
    {
      std::string tmp = "FixELimit";
      const int max_type = force->registry.max_type();
      MatrixProperty *coeffRestProp = force->registry.getMatrixProperty(COEFFICIENT_RESTITUTION, tmp.c_str());
      double **coeffRest = coeffRestProp->data;
      for (int i = 1; i < max_type + 1; i++)
      {
        for (int j = 1; j < max_type + 1; j++)
        {
          coravg += coeffRest[i][j];
          id++;
        }
      }
    }
    coravg /= (id > 0) ? id : 1.0;

    std::string additional = "";

    if (fix_ms)
    {
      int nbodyall = fix_ms->n_body_all();
      additional.append(" nbody: ");
      additional.append(std::to_string(nbodyall));
      additional.append(" ");
    }

    if (comm->me == 0)
    {
      std::cout << std::scientific << std::setprecision(3) << "pre-run, step: " << update->ntimestep << ", ke/mass=" << cke / global_sum_mass
                << ", natoms=" << atom->natoms << ", overlap=" << global_minOverlap << ", cor=" << coravg << ", dt=" << update->dt << additional
                << "\n"; //<<", cint="<<cke_icheck <<"\n";
      // std::cout << "ckem: ";
      // for (int i=0; i<10; i++) std::cout <<cke_mass[i]<<", ";
      // std::cout << "\n";
    }
    if (cke / global_sum_mass < threshold)
      check_exit = 1;
  }
  // }
  int resvv = 0;
  // if (update->ntimestep % 10000 == 0) cke_icheck++;
  // // shift memory and record new one
  // if (update->ntimestep % (10)==0) {
  //   double cke_mass_tmp[10];
  //   for (int i=0; i<10; i++) cke_mass_tmp[i] = cke_mass[i];
  //   cke_mass[0] = cke/global_sum_mass;
  //   for (int i=1; i<10; i++) cke_mass[i] = cke_mass_tmp[i-1];

  // if (update->ntimestep>=1000) {
  //   // compute derivative d Ke/m / dt, if negative add 1
  //   double dtf1 = (cke_mass[0] - cke_mass[1]);
  //   double dtf2 = (cke_mass[0] - 4./3. * cke_mass[1] + 1./3.*cke_mass[2])*1.5;
  //   double dtf3 = (cke_mass[0] - 18./11. * cke_mass[1] + 9./11.*cke_mass[2] - 2./11.*cke_mass[3])*11./6.;
  //   // if (comm->me == 0 && update->ntimestep % 1000 == 0) std::cout << "derivatives: "<< dtf1<<", "<<dtf2<<", "<<dtf3<< "\n";
  //   // if ((dtf1<0. && dtf2<0. && dtf3<0.) &&
  //   //     (abs(dtf1)<0. && abs(dtf2)<0. && abs(dtf3)<1e-1)
  //   //     )
  //   if (dtf1<0. && dtf2<0. && dtf3<0.)
  //     check_n_int++;
  //   if (dtf1>0. && dtf2>0. && dtf3>0.)
  //     check_n_int--;
  //   if (check_n_int>10) {
  //     check_n_int = 0;
  //     cke_icheck++;
  //   }
  //   check_n_int = std::max(check_n_int,0);
  //   cke_icheck = std::max(std::min(cke_icheck,check_int),2);
  // }
  // }

  // if (// update->ntimestep<1000 &&
  //     update->ntimestep % (cke_icheck)==0) resvv=1;
  // // afterwards, do analysis
  // // for (int i=0; i<10; i++)
  // //   if (update->ntimestep>=i*1000 && update->ntimestep<(i+1)*1000 && update->ntimestep % (i+2)==0)  resvv=1;
  if (update->ntimestep % (check_int) == 0)
    resvv = 1;

  // if (update->ntimestep % check_int == 0) {
  if (resvv)
  {
    reset_vel();
    if (check_exit)
    {
      {
        // reset variable parameters
        update->dt = end_dt;

        std::string tmp = "FixELimit";
        const int max_type = force->registry.max_type();
        MatrixProperty *coeffRestProp = force->registry.getMatrixProperty(COEFFICIENT_RESTITUTION, tmp.c_str());
        double **coeffRest = coeffRestProp->data;
        int id = 0;
        for (int i = 1; i < max_type + 1; i++)
        {
          for (int j = 1; j < max_type + 1; j++)
          {
            coeffRest[i][j] = cor0[id++];
          }
        }
      }
      SignalHandler::request_exit = true;
    }
  }

  double **v = atom->v;
  double ave_velocity[3] = {0.0};
  double total_mass = 0.0;

  if (!fix_ms)
  {

    for (int i = 0; i < atom->nlocal; i++)
    {
      ave_velocity[0] += atom->rmass[i] * v[i][0];
      ave_velocity[1] += atom->rmass[i] * v[i][1];
      ave_velocity[2] += atom->rmass[i] * v[i][2];

      total_mass += atom->rmass[i];
    }

    double global_ave_velocity[3];
    MPI_Allreduce(&ave_velocity[0], &global_ave_velocity[0],
                  3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double global_total_mass;
    MPI_Allreduce(&total_mass, &global_total_mass,
                  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (int i = 0; i < 3; i++)
    {
      global_ave_velocity[i] = global_ave_velocity[i] / global_total_mass;
    }

    for (int i = 0; i < atom->nlocal; i++)
    {
      v[i][0] -= global_ave_velocity[0];
      v[i][1] -= global_ave_velocity[1];
      v[i][2] -= global_ave_velocity[2];
    }
  }
  else if (fix_ms)
  {

    ave_velocity[0] = 0.0;
    ave_velocity[1] = 0.0;
    ave_velocity[2] = 0.0;
    total_mass = 0.0;

    for (int i = 0; i < atom->nlocal; i++)
    {
      if (fix_ms->belongs_to(i) == -1)
      {
        ave_velocity[0] += atom->rmass[i] * v[i][0];
        ave_velocity[1] += atom->rmass[i] * v[i][1];
        ave_velocity[2] += atom->rmass[i] * v[i][2];

        total_mass += atom->rmass[i];
      }
    }

    double vcm[3];
    for (int i = 0; i < fix_ms->n_body(); i++)
    {

      fix_ms->data().vcm(vcm, i);
      ave_velocity[0] += fix_ms->data().mass(i) * vcm[0];
      ave_velocity[1] += fix_ms->data().mass(i) * vcm[1];
      ave_velocity[2] += fix_ms->data().mass(i) * vcm[2];

      total_mass += fix_ms->data().mass(i);
    }

    double global_ave_velocity[3];
    MPI_Allreduce(&ave_velocity[0], &global_ave_velocity[0],
                  3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double global_total_mass;
    MPI_Allreduce(&total_mass, &global_total_mass,
                  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (int i = 0; i < 3; i++)
    {
      global_ave_velocity[i] = global_ave_velocity[i] / global_total_mass;
    }

    for (int i = 0; i < atom->nlocal; i++)
    {
      if (fix_ms->belongs_to(i) == -1)
      {
        v[i][0] -= global_ave_velocity[0];
        v[i][1] -= global_ave_velocity[1];
        v[i][2] -= global_ave_velocity[2];
      }
    }

    for (int i = 0; i < fix_ms->n_body(); i++)
    {

      fix_ms->data().vcm(vcm, i);
      vcm[0] -= global_ave_velocity[0];
      vcm[1] -= global_ave_velocity[1]; // + 0.2 * domain->yprd * ssr;
      vcm[2] -= global_ave_velocity[2];
      fix_ms->data().set_v_body(i, vcm);
    }
    fix_ms->set_xv();
  }
}
