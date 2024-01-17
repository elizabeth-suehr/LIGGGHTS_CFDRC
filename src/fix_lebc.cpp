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
------------------------------------------------------------------------- */

/*------------------------------------------------------------------------
    This is code for the Lees Edwards Boundary Condition for Multisphere
    and single sphere, written by Elizabeth Suehr
    elizabeth.suehr@gmail.com
    emsuehr@ucdavis.edu //Probably gone after 2025-2026

    This file is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with this code.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    I do not know how licensing work, please email me if you see anything
    wrong so I do not get in trouble :)
------------------------------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include "fix_lebc.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "pair.h"

#include <random>
#include <fstream>
#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace std;

/* ---------------------------------------------------------------------- */

FixLEBC::FixLEBC(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg < 9)
    error->fix_error(FLERR, this, "Not enough arguments, include a float value for the shear strain rate, and then true or false for if the output is dimensional\n Example:fix     lebounary all lebc 100.0 true gtemp 0.01 ave_reset 50000 body_data curl1\n");

  ssr = force->numeric(FLERR, arg[3]);
  domain->ssr = ssr;
  // domain->lebc_displacement = domain->xprd * 0.0001;

  if (strcmp(arg[4], "false") == 0)
  {
    isDimensional = false;
  }

  if (strcmp(arg[5], "gtemp") == 0)
  {
    gtemp_distribution = force->numeric(FLERR, arg[6]);
  }

  if (strcmp(arg[7], "ave_reset") == 0)
  {
    ave_count_reset = force->inumeric(FLERR, arg[8]);
  }

  if (strcmp(arg[9], "body_data") == 0)
  {
    body_data_name = arg[10];
    ave_count_reset = force->inumeric(FLERR, arg[11]);
  }

  if (domain->nonperiodic != 0 ||
      domain->xperiodic != 1 ||
      domain->yperiodic != 1 ||
      domain->zperiodic != 1)
  {
    error->fix_error(FLERR, this, "Set domain to be periodic in all directions \n");
  }

  if (!comm->ghost_velocity)
  {
    error->fix_error(FLERR, this, "ghost particles should transfer velocities \n");
  }

  if (domain->triclinic == 1)
  {
    error->fix_error(FLERR, this, "Don't use triclinic \n");
  }

  fix_ms = static_cast<FixMultisphere *>(modify->find_fix_style("multisphere", 0));
  force->pair->pre_init_lebc();

  atom->add_callback(2);
  comm_border = 1;

  if (fix_ms)
  {
    init_total_bodies = fix_ms->n_body_all();
  }
}
/* ---------------------------------------------------------------------- */

FixLEBC::~FixLEBC()
{
  atom->delete_callback(id, 2);
}

/* ---------------------------------------------------------------------- */

int FixLEBC::setmask()
{
  int mask = 0;
  // mask |= INITIAL_INTEGRATE;
  // mask |= PRE_EXCHANGE;
  mask |= PRE_NEIGHBOR;
  // mask |= PRE_FORCE;
  // mask |= PRE_FINAL_INTEGRATE;
  // mask |= FINAL_INTEGRATE;
  mask |= POST_INTEGRATE;
  mask |= POST_FORCE;
  // mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLEBC::init()
{
  domain->lebc = 1;

  double **x = atom->x;
  double **v = atom->v;

  // random device class instance, source of 'true' randomness for initializing random seed
  random_device rd;

  // Mersenne twister PRNG, initialized with seed from previous random device instance
  // mt19937 gen(rd());
  mt19937 gen(10);

  double ave_velocity[3] = {0.0};
  double total_mass = 0.0;

  if (!fix_ms)
  {
    for (int i = 0; i < atom->nlocal; i++)
    {

      double velocity = (x[i][1] - 0.5 * domain->yprd) * ssr;

      // instance of class std::normal_distribution with specific mean and stddev
      std::normal_distribution<float> d_x(velocity, gtemp_distribution * domain->yprd * ssr);
      std::normal_distribution<float> d_y(0.0, gtemp_distribution * domain->yprd * ssr);
      std::normal_distribution<float> d_z(0.0, gtemp_distribution * domain->yprd * ssr);

      v[i][0] = d_x(gen);
      v[i][1] = d_y(gen);
      v[i][2] = d_z(gen);
    }

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

    for (int i = 0; i < atom->nlocal; i++)
    {
      if (fix_ms->belongs_to(i) == -1)
      {
        double velocity = (x[i][1] - 0.5 * domain->yprd) * ssr;

        // instance of class std::normal_distribution with specific mean and stddev
        std::normal_distribution<float> d_x(velocity, gtemp_distribution * domain->yprd * ssr);
        std::normal_distribution<float> d_y(0.0, gtemp_distribution * domain->yprd * ssr);
        std::normal_distribution<float> d_z(0.0, gtemp_distribution * domain->yprd * ssr);

        v[i][0] = d_x(gen);
        v[i][1] = d_y(gen);
        v[i][2] = d_z(gen);
      }
    }

    for (int ibody = 0; ibody < fix_ms->n_body(); ibody++)
    {

      double xcm[3];

      fix_ms->data().xcm(xcm, ibody);
      double vel = (xcm[1] - 0.5 * domain->yprd) * ssr;

      // instance of class std::normal_distribution with specific mean and stddev
      std::normal_distribution<float> d_x(vel, gtemp_distribution * domain->yprd * ssr);
      std::normal_distribution<float> d_y(0.0, gtemp_distribution * domain->yprd * ssr);
      std::normal_distribution<float> d_z(0.0, gtemp_distribution * domain->yprd * ssr);

      double velocity[3] = {d_x(gen), d_y(gen), d_z(gen)};

      fix_ms->data().set_v_body(ibody, velocity);
    }

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
  // TODO might not have the averaging of mixed single sphere and multisphere done correctly
}

// NOT CALLED, used for debugging
void FixLEBC::initial_integrate(int vflag)
{
  int xbox, ybox, zbox;
  tagint *image = atom->image;

  for (int i = 0; i < atom->nlocal + atom->nghost; i++)
  {

    xbox = (image[i] & IMGMASK) - IMGMAX;
    ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (image[i] >> IMG2BITS) - IMGMAX;

    if (xbox > 1 || ybox > 1 || xbox < -1 || ybox < -1)
    {
      if (i >= atom->nlocal)
      {
        std::cout << "lebcg " << xbox << " " << ybox << "\n";
      }
      else
      {
        std::cout << "lebcr " << xbox << " " << ybox << "\n";
      }
    }
  }
}

int FixLEBC::pack_border(int n, int *list, double *buf, int *pbc)
{

  tagint *image = atom->image;
  int i, j, m = 0;
  for (i = 0; i < n; i++)
  {
    j = list[i];

    buf[m++] = static_cast<double>(pbc[1]);
    // buf[m++] = static_cast<double>(image[j]);
  }
  return m;
}

int FixLEBC::unpack_border(int n, int first, double *buf)
{
  int i, m, last;

  double **x = atom->x;
  double **v = atom->v;

  tagint *image = atom->image;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
  {
    int pbc_1 = static_cast<int>(buf[m++]);

    x[i][0] += pbc_1 * domain->lebc_displacement;
    v[i][0] += pbc_1 * domain->yprd * domain->ssr;

    if (pbc_1)
    {

      if (x[i][0] > domain->xprd)
      {
        x[i][0] -= domain->xprd;
      }
      if (x[i][0] < 0.0)
      {
        x[i][0] += domain->xprd;
      }
    }

    // int im = static_cast<int>(buf[m++]);
    // image[i] = im;
  }

  return m;
}

void FixLEBC::pre_neighbor()
{
  // used to test comm setup without multispheres
  // very useful!!!

  // if (!fix_ms)
  // {
  //   forward_comm();
  //   reverse_comm();
  // }
}

void FixLEBC::forward_comm()
{
  comm->forward_comm_fix(this);
}

int FixLEBC::pack_comm_f_test(int n, int *list, double *buf, int pbc_flag, int *pbc)
{

  int i, j, m = 0;
  for (i = 0; i < n; i++)
  {
    j = list[i];

    buf[m++] = static_cast<double>(pbc[1]);
  }
  return 1;
}

int FixLEBC::pack_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{

  return pack_comm_f_test(n, list, buf, pbc_flag, pbc);
}

void FixLEBC::unpack_comm_f_test(int n, int first, double *buf)
{
  int i, m, last;

  double **x = atom->x;
  double **v = atom->v;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
  {
    int pbc_1 = static_cast<int>(buf[m++]);

    if (pbc_1)
    {
    }
  }
}

void FixLEBC::unpack_comm(int n, int first, double *buf)
{
  unpack_comm_f_test(n, first, buf);
  return;
}

void FixLEBC::reverse_comm()
{
  comm->reverse_comm_fix(this);
}

/* ----------------------------------------------------------------------
   pack reverse comm
------------------------------------------------------------------------- */

int FixLEBC::pack_reverse_comm(int n, int first, double *buf)
{

  return pack_reverse_comm_f_test(n, first, buf);
}

/* ---------------------------------------------------------------------- */

int FixLEBC::pack_reverse_comm_f_test(int n, int first, double *buf)
{
  int i, m, last, tag, flag;

  double **x = atom->x;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
  {
    buf[m++] = static_cast<double>(0);
  }
  return 1; // one is length of buff per atom exchanged
}

/* ---------------------------------------------------------------------- */

void FixLEBC::setup(int vflag)
{
}

/* ---------------------------------------------------------------------- */

void FixLEBC::post_integrate()
{
  domain->update_lebc_displacement();
}

/* ---------------------------------------------------------------------- */

void FixLEBC::post_force(int vflag)
{

  double domain_volume = domain->xprd * domain->yprd * domain->zprd;

  double **x = atom->x;
  double **v = atom->v;

  double ave_velocity[3] = {0.0};

  double total_mass = 0.0;

  for (int i = 0; i < atom->nlocal; i++)
  {
    if (!fix_ms)
    {
      ave_velocity[0] += atom->rmass[i] * v[i][0];
      ave_velocity[1] += atom->rmass[i] * v[i][1];
      ave_velocity[2] += atom->rmass[i] * v[i][2];

      total_mass += atom->rmass[i];
    }
    else if (fix_ms->belongs_to(i) == -1)
    {
      ave_velocity[0] += atom->rmass[i] * v[i][0];
      ave_velocity[1] += atom->rmass[i] * v[i][1];
      ave_velocity[2] += atom->rmass[i] * v[i][2];

      total_mass += atom->rmass[i];
    }
  }

  if (fix_ms)
  {
    double vcm[3];
    for (int i = 0; i < fix_ms->n_body(); i++)
    {

      fix_ms->data().vcm(vcm, i);
      ave_velocity[0] += fix_ms->data().mass(i) * vcm[0];
      ave_velocity[1] += fix_ms->data().mass(i) * vcm[1];
      ave_velocity[2] += fix_ms->data().mass(i) * vcm[2];

      total_mass += fix_ms->data().mass(i);
    }
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
    // global_ave_velocity[i] = 0.0;
  }

  double local_kst[6] = {0};

  // TODO? Right now, we are subtracting ave_velocity, even though we shouldn't need to because
  // we subtract it at the begining, maybe removing in future? Side note, The atoms are still drifting up and down, which they shouldn't
  // so maybe we need to keep it
  for (int i = 0; i < atom->nlocal; i++)
  {
    if (!fix_ms)
    {
      double local_velocity = (x[i][1] - 0.5 * domain->yprd) * ssr;
      local_kst[0] += atom->rmass[i] * (v[i][0] - global_ave_velocity[0] - local_velocity) * (v[i][0] - global_ave_velocity[0] - local_velocity); // xx
      local_kst[1] += atom->rmass[i] * (v[i][1] - global_ave_velocity[1]) * (v[i][1] - global_ave_velocity[1]);                                   // yy
      local_kst[2] += atom->rmass[i] * (v[i][2] - global_ave_velocity[2]) * (v[i][2] - global_ave_velocity[2]);                                   // zz
      local_kst[3] += atom->rmass[i] * (v[i][0] - global_ave_velocity[0] - local_velocity) * (v[i][1] - global_ave_velocity[1]);                  // xy
      local_kst[4] += atom->rmass[i] * (v[i][0] - global_ave_velocity[0] - local_velocity) * (v[i][2] - global_ave_velocity[2]);                  // xz
      local_kst[5] += atom->rmass[i] * (v[i][1] - global_ave_velocity[1]) * (v[i][2] - global_ave_velocity[2]);                                   // yz
    }
    else if (fix_ms->belongs_to(i) == -1)
    {
      double local_velocity = (x[i][1] - 0.5 * domain->yprd) * ssr;
      local_kst[0] += atom->rmass[i] * (v[i][0] - global_ave_velocity[0] - local_velocity) * (v[i][0] - global_ave_velocity[0] - local_velocity); // xx
      local_kst[1] += atom->rmass[i] * (v[i][1] - global_ave_velocity[1]) * (v[i][1] - global_ave_velocity[1]);                                   // yy
      local_kst[2] += atom->rmass[i] * (v[i][2] - global_ave_velocity[2]) * (v[i][2] - global_ave_velocity[2]);                                   // zz
      local_kst[3] += atom->rmass[i] * (v[i][0] - global_ave_velocity[0] - local_velocity) * (v[i][1] - global_ave_velocity[1]);                  // xy
      local_kst[4] += atom->rmass[i] * (v[i][0] - global_ave_velocity[0] - local_velocity) * (v[i][2] - global_ave_velocity[2]);                  // xz
      local_kst[5] += atom->rmass[i] * (v[i][1] - global_ave_velocity[1]) * (v[i][2] - global_ave_velocity[2]);                                   // yz
    }
  }

  if (fix_ms)
  {

    double xcm[3];
    double vcm[3];
    for (int i = 0; i < fix_ms->n_body(); i++)
    {

      fix_ms->data().xcm(xcm, i);
      fix_ms->data().vcm(vcm, i);
      double local_velocity = (xcm[1] - 0.5 * domain->yprd) * ssr;

      local_kst[0] += fix_ms->data().mass(i) * (vcm[0] - global_ave_velocity[0] - local_velocity) * (vcm[0] - global_ave_velocity[0] - local_velocity); // xx
      local_kst[1] += fix_ms->data().mass(i) * (vcm[1] - global_ave_velocity[1]) * (vcm[1] - global_ave_velocity[1]);                                   // yy
      local_kst[2] += fix_ms->data().mass(i) * (vcm[2] - global_ave_velocity[2]) * (vcm[2] - global_ave_velocity[2]);                                   // zz
      local_kst[3] += fix_ms->data().mass(i) * (vcm[0] - global_ave_velocity[0] - local_velocity) * (vcm[1] - global_ave_velocity[1]);                  // xy
      local_kst[4] += fix_ms->data().mass(i) * (vcm[0] - global_ave_velocity[0] - local_velocity) * (vcm[2] - global_ave_velocity[2]);                  // xz
      local_kst[5] += fix_ms->data().mass(i) * (vcm[1] - global_ave_velocity[1]) * (vcm[2] - global_ave_velocity[2]);                                   // yz
    }
  }

  double global_kst[6] = {0.0};

  MPI_Allreduce(&local_kst[0], &global_kst[0],
                6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (int i = 0; i < 6; i++)
  {
    global_kst[i] = global_kst[i] / domain_volume;
  }

  for (int i = 0; i < 6; i++)
  {
    kinetic_stress_tensor[i] = (kinetic_stress_tensor[i] * ave_count + global_kst[i]) / ((float)(ave_count + 1));
  }

  double local_cst[6] = {0.0};
  if (force->pair)
  {
    double **force_length_atom = force->pair->force_length;

    for (int i = 0; i < atom->nlocal; i++)
    {
      for (int j = 0; j < 6; j++)
      {
        local_cst[j] += force_length_atom[i][j];
      }
    }
  }

  if (fix_ms)
  {
    for (int i = 0; i < atom->nlocal + atom->nghost; i++)
    {
      for (int j = 0; j < 6; j++)
      {
        local_cst[j] += fix_ms->force_length[i][j];
      }
    }
  }

  double global_cst[6] = {0.0};
  MPI_Allreduce(&local_cst[0], &global_cst[0], 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (int i = 0; i < 6; i++)
  {
    global_cst[i] = global_cst[i] / domain_volume;
  }

  for (int i = 0; i < 6; i++)
  {
    collision_stress_tensor[i] = (collision_stress_tensor[i] * ave_count + global_cst[i]) / ((float)(ave_count + 1));
  }

  ave_count++;

  if (ave_count % 2500 == 0)
  {
    double equiv_vol_diameter = 1.0;
    if (!isDimensional)
    {
      if (fix_ms)
      {
        equiv_vol_diameter = cbrt(6.0 * fix_ms->data().volume(0) / M_PI);
      }
      else
      {
        equiv_vol_diameter = 2.0 * atom->radius[0];
      }
    }

    double scale = 1.0;
    if (!isDimensional)
    {
      scale = 1.0 / (2500.0 * 100.0 * 100.0 * (equiv_vol_diameter) * (equiv_vol_diameter));
      ;
    }

    if (comm->me == 0)
    {
      if (fix_ms)
        cout << update->ntimestep << " " << ave_count << " body_count: " << fix_ms->n_body_all() << " atom_count: " << atom->natoms << " " << kinetic_stress_tensor[0] * scale << " " << kinetic_stress_tensor[1] * scale << " " << kinetic_stress_tensor[3] * scale << " " << collision_stress_tensor[0] * scale << " " << collision_stress_tensor[1] * scale << " " << collision_stress_tensor[3] * scale << "\n";
      else
        cout << update->ntimestep << " " << ave_count << " atom_count: " << atom->natoms << " " << kinetic_stress_tensor[0] * scale << " " << kinetic_stress_tensor[1] * scale << " " << kinetic_stress_tensor[3] * scale << " " << collision_stress_tensor[0] * scale << " " << collision_stress_tensor[1] * scale << " " << collision_stress_tensor[3] * scale << "\n";
    }
  }

  if (ave_count == ave_count_reset)
  {

    if (comm->me == 0)
    {
      double equiv_vol_diameter = 1.0;
      if (!isDimensional)
      {
        if (fix_ms)
        {
          equiv_vol_diameter = cbrt(6.0 * fix_ms->data().volume(0) / M_PI);
        }
        else
        {
          equiv_vol_diameter = 2.0 * atom->radius[0];
        }
      }

      double scale = 1.0;
      if (!isDimensional)
      {
        scale = 1.0 / (2500.0 * 100.0 * 100.0 * (equiv_vol_diameter) * (equiv_vol_diameter));
      }
      cout << update->ntimestep << " " << kinetic_stress_tensor[0] * scale << " " << kinetic_stress_tensor[1] * scale << " " << kinetic_stress_tensor[3] * scale << " " << collision_stress_tensor[0] * scale << " " << collision_stress_tensor[1] * scale << " " << collision_stress_tensor[3] * scale << "\n";
      fprintf(logfile, "stress: %f %f %f %f %f %f\n", kinetic_stress_tensor[0] * scale, kinetic_stress_tensor[1] * scale, kinetic_stress_tensor[3] * scale, collision_stress_tensor[0] * scale, collision_stress_tensor[1] * scale, collision_stress_tensor[3] * scale);
      fflush(logfile);
    }
    ave_count = 0;
    for (int i = 0; i < 6; i++)
    {
      kinetic_stress_tensor[i] = 0.0;
      collision_stress_tensor[i] = 0.0;
    }
  }

  // if (fix_ms)
  // {
  //   int nbody = fix_ms->n_body();
  //   for (int ibody = 0; ibody < nbody; ibody++)
  //   {
  //     // int image = fix_ms->data().imagebody_.get(ibody);
  //     int remap[4];
  //     fix_ms->data().remap_values(remap,ibody);
  //     // double vcm[3];
  //     // fix_ms->data().vcm(vcm, ibody);
  //     // double xcm[3];
  //     // fix_ms->data().xcm(xcm, ibody);
  //     // std::cout << vcm[0] << " " << vcm[1] << " " << vcm[2] << "\n";
  //     // std::cout << xcm[0] << " " << xcm[1] << " " << xcm[2] << "\n";

  //     // std::cout << remap[0] << " " << remap[1] << " " << remap[2] << " " << remap[3] << "\n";
  //   }
  // }
  // fix_ms->check_bodies_atoms("Fixlebc::post_force()");

  // if (fix_ms)
  // {
  //   if (fix_ms->n_body_all() != init_total_bodies)
  //   {
  //     std::cout << "n_body_all changed to: " << fix_ms->n_body_all() << " from: " << init_total_bodies << "\n";
  //     fprintf(logfile, "n_body_all changed to: %i from: %i",fix_ms->n_body_all(),  init_total_bodies);
  //     fflush(logfile);

  //     MPI_Abort(0,42);
  //   }
  // }
  if (save_count == save_count_reset)
  {
    print_body_data();
    save_count = 0;
  }
  save_count++;
}

void FixLEBC::print_body_data()
{
  std::vector<double> body_rotations_total(fix_ms->n_body_all() * 4);
  std::vector<double> body_rotations_local(fix_ms->n_body() * 4);

  std::vector<double> body_positions_total(fix_ms->n_body_all() * 3);
  std::vector<double> body_positions_local(fix_ms->n_body() * 3);

  std::vector<int> body_tag_total(fix_ms->n_body_all());
  std::vector<int> body_tag_local(fix_ms->n_body());

  for (int ibody = 0; ibody < fix_ms->n_body(); ibody++)
  {

    double quat[4];
    double pos[3];

    fix_ms->data().quat(quat, ibody);
    fix_ms->data().xcm(pos, ibody);

    body_rotations_local.push_back(quat[0]);
    body_rotations_local.push_back(quat[1]);
    body_rotations_local.push_back(quat[2]);
    body_rotations_local.push_back(quat[3]);

    body_positions_local.push_back(pos[0]);
    body_positions_local.push_back(pos[1]);
    body_positions_local.push_back(pos[2]);

    body_tag_local.push_back(fix_ms->data().tag(ibody));
  }

  MPI_Allgather(&body_rotations_local, fix_ms->n_body() * 4, MPI_DOUBLE, &body_rotations_total, fix_ms->n_body_all() * 4, MPI_DOUBLE, world);
  MPI_Allgather(&body_positions_local, fix_ms->n_body() * 3, MPI_DOUBLE, &body_positions_total, fix_ms->n_body_all() * 3, MPI_DOUBLE, world);
  MPI_Allgather(&body_tag_local, fix_ms->n_body(), MPI_INT, &body_tag_total, fix_ms->n_body_all(), MPI_INT, world);

  if (comm->me == 0)
  {
    ofstream myfile;
    string filename = body_data_name + "/" + std::to_string(update->ntimestep) + ".txt";

    myfile.open(filename);

    for (int ibody = 0; ibody < fix_ms->n_body_all(); ibody++)
    {
      myfile << ibody << " " << body_tag_total[ibody] << " "
             << body_rotations_total[ibody * 4 + 0] << " "
             << body_rotations_total[ibody * 4 + 1] << " "
             << body_rotations_total[ibody * 4 + 2] << " "
             << body_rotations_total[ibody * 4 + 3]
             << " " << body_positions_local[ibody * 3 + 0]
             << " " << body_positions_local[ibody * 3 + 1]
             << " " << body_positions_local[ibody * 3 + 2] << "\n";
    }

    myfile.close();
  }
}
