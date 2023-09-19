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
/*------------------------------------------------------------------------
    This is code has been edited for the Lees Edwards Boundary Condition for
    Multisphere and single sphere, edited by Elizabeth Suehr
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

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fix_multisphere.h"
#include "domain_wedge.h"
#include "math_extra.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "respa.h"
#include "modify.h"
#include "group.h"
#include "comm.h"
#include "force.h"
#include "output.h"
#include "memory.h"
#include "error.h"
#include "fix_property_atom.h"
#include "fix_template_multisphere.h"
#include "neighbor.h"
#include "fix_gravity.h"
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "atom_vec.h"
#include "math_extra_liggghts.h"
#include "math_const.h"


#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define TOLERANCE 1.0e-6
#define EPSILON 1.0e-7
#define MAXJACOBI 50
#define DELTA_GROW 10000

enum
{
    LOOP_LOCAL,
    LOOP_ALL
};

#ifndef LMP_MULTISPHERE_PARALLEL_FLAG_
typedef Multisphere MultisphereParallel;
#endif // LMP_MULTISPHERE_PARALLEL_FLAG_

/* ---------------------------------------------------------------------- */

FixMultisphere::FixMultisphere(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg),
                                                                    do_modify_body_forces_torques_(false),
                                                                    multisphere_(*(new MultisphereParallel(lmp))),
                                                                    fix_corner_ghost_(0),
                                                                    fix_delflag_(0),
                                                                    fix_existflag_(0),
                                                                    fix_volumeweight_ms_(0),
                                                                    use_volumeweight_ms_(true),
                                                                    fix_gravity_(0),
                                                                    fw_comm_flag_(MS_COMM_UNDEFINED),
                                                                    rev_comm_flag_(MS_COMM_UNDEFINED),
                                                                    body_(NULL),
                                                                    force_length(NULL),
                                                                    displace_(NULL),
                                                                    ntypes_(0),
                                                                    Vclump_(0),
                                                                    allow_group_and_set_(false),
                                                                    allow_heatsource_(false),
                                                                    CAdd_(0.),
                                                                    fluidDensity_(0.),
                                                                    concave_(false),
                                                                    add_dragforce_(true)
{

    if (0 == strcmp(style, "concave"))
    {
        concave_ = true;
        int strln = strlen("multisphere") + 1;
        delete[] style;
        style = new char[strln];
        strcpy(style, "multisphere");
    }

    int iarg = 3;

    bool hasargs = true;
    while (iarg < narg && hasargs)
    {
        hasargs = false;
        if (strcmp(arg[iarg], "allow_group_and_set") == 0)
        {
            if (narg < iarg + 2)
                ms_error(FLERR, "not enough arguments for 'allow_group_and_set'");
            if (0 == strcmp(arg[iarg + 1], "yes"))
                allow_group_and_set_ = true;
            else if (0 == strcmp(arg[iarg + 1], "no"))
                allow_group_and_set_ = false;
            else
                ms_error(FLERR, "expecting 'yes' or 'no' after 'allow_group_and_set'");
            iarg += 2;
            hasargs = true;
        }
        else if (strcmp(arg[iarg], "allow_heatsource") == 0)
        {
            if (narg < iarg + 2)
                ms_error(FLERR, "not enough arguments for 'allow_heatsource'");
            if (strcmp(arg[iarg + 1], "yes") == 0)
                allow_heatsource_ = true;
            else if (strcmp(arg[iarg + 1], "no"))
                allow_heatsource_ = false;
            else
                ms_error(FLERR, "expecting 'yes' or 'no' after 'allow_heatsource'");
            iarg += 2;
            hasargs = true;
        }
        else if (strcmp(arg[iarg], "CAddRhoFluid") == 0)
        {
            if (narg < iarg + 3)
                ms_error(FLERR, "not enough arguments for 'CAddRhoFluid'. You must specify the added mass coefficient AND the fluid density");
            CAdd_ = atof(arg[iarg + 1]);
            fluidDensity_ = atof(arg[iarg + 2]);
            fprintf(screen, "cfd_coupling_force_ms_implicit will consider added mass with CAdd = %g, fluidDensity: %g\n",
                    CAdd_, fluidDensity_);
            iarg += 3;
            hasargs = true;
        }
        else if (0 == strcmp(style, "multisphere") || 0 == strcmp(style, "multisphere/advanced"))
        {
            char *errmsg = new char[strlen(arg[iarg]) + 50];
            sprintf(errmsg, "unknown keyword or wrong keyword order: %s", arg[iarg]);
            ms_error(FLERR, errmsg);
            delete[] errmsg;
        }
    }

    if (atom->molecular == 1)
        ms_error(FLERR, "Must NOT use a hybrid sphere/molecular atom style with fix multisphere (use sphere only)");

    atom->molecule_flag = 1;
    grow_arrays(atom->nmax);

    char **modarg;
    modarg = new char *[3];
    modarg[2] = new char[50];
    modarg[0] = (char *)"exclude";
    modarg[1] = (char *)"molecule";
    strcpy(modarg[2], arg[1]);
    neighbor->modify_params(3, modarg);
    delete[] modarg[2];
    delete[] modarg;

    restart_global = 1;
    restart_peratom = 1;
    restart_pbc = 1;
    atom->add_callback(0);
    atom->add_callback(1);

    if (accepts_restart_data_from_style)
        delete[] accepts_restart_data_from_style;
    accepts_restart_data_from_style = new char[21];
    sprintf(accepts_restart_data_from_style, "multisphere/advanced");

    // fix handles properties that need to be initialized at particle creation
    create_attribute = 1;

    force_reneighbor = 1;
    next_reneighbor = -1;

    // is now local data, not global
    local_flag = 1;

    size_local_rows = 0;
    size_local_cols = 12; // 0 = vector, N = columns in local array
    local_freq = 1;

    size_peratom_cols = 0;

    vector_flag = 1;
    size_vector = 0; // no bodies present at creation

    global_freq = 1;
    extarray = 0;

    comm_forward = 7;

    comm_reverse = 10;

    if (atom->quaternion)
        comm_reverse += 4;

    lebc_highest_body_id = 0;
}

/* ---------------------------------------------------------------------- */

FixMultisphere::~FixMultisphere()
{
    atom->delete_callback(id, 0);
    atom->delete_callback(id, 1);

    delete &multisphere_;

    memory->destroy(displace_);

    if (accepts_restart_data_from_style)
    {
        delete[] accepts_restart_data_from_style;
        accepts_restart_data_from_style = 0;
    }
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::ms_error(const char *file, int line, char const *errmsg)
{
    if (concave_)
        error->fix_error(file, line, this, "concave", errmsg);
    else
        error->fix_error(file, line, this, errmsg);
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::post_create()
{

    if (!fix_corner_ghost_)
    {
        const char *fixarg[9];
        fixarg[0] = "cornerghost";
        fixarg[1] = "all";
        fixarg[2] = "property/atom";
        fixarg[3] = "cornerghost";
        fixarg[4] = "scalar";
        fixarg[5] = "no"; // restart
        fixarg[6] = "no"; // communicate ghost
        fixarg[7] = "no"; // communicate rev
        fixarg[8] = "0.";
        fix_corner_ghost_ = modify->add_fix_property_atom(9, const_cast<char **>(fixarg), style);
    }

    if (!fix_delflag_)
    {
        const char *fixarg[9];
        fixarg[0] = "delflag";
        fixarg[1] = "all";
        fixarg[2] = "property/atom";
        fixarg[3] = "delflag";
        fixarg[4] = "scalar";
        fixarg[5] = "yes"; // restart
        fixarg[6] = "no";  // communicate ghost
        fixarg[7] = "yes"; // communicate rev
        fixarg[8] = "0.";
        fix_delflag_ = modify->add_fix_property_atom(9, const_cast<char **>(fixarg), style);
    }

    if (!fix_existflag_)
    {
        const char *fixarg[9];
        fixarg[0] = "existflag";
        fixarg[1] = "all";
        fixarg[2] = "property/atom";
        fixarg[3] = "existflag";
        fixarg[4] = "scalar";
        fixarg[5] = "no";  // restart
        fixarg[6] = "no";  // communicate ghost
        fixarg[7] = "yes"; // communicate rev
        fixarg[8] = "1.";
        fix_existflag_ = modify->add_fix_property_atom(9, const_cast<char **>(fixarg), style);
    }

    if (!fix_volumeweight_ms_ && use_volumeweight_ms_)
    {
        const char *fixarg[9];
        fixarg[0] = "volumeweight_ms";
        fixarg[1] = "all";
        fixarg[2] = "property/atom";
        fixarg[3] = "volumeweight_ms";
        fixarg[4] = "scalar";
        fixarg[5] = "yes"; // restart
        fixarg[6] = "yes"; // communicate ghost
        fixarg[7] = "no";  // communicate rev
        fixarg[8] = "1.";
        fix_volumeweight_ms_ = modify->add_fix_property_atom(9, const_cast<char **>(fixarg), style);
    }

    if (modify->have_restart_data(this))
    {
        evflag = 0;
        set_xv(LOOP_LOCAL);
    }
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::pre_delete(bool unfixflag)
{
    if (unfixflag)
        ms_error(FLERR, "this fix may not be unfixed as it holds "
                        "all the internal data for multi-spheres");
}

/* ---------------------------------------------------------------------- */

int FixMultisphere::setmask()
{
    int mask = 0;
    mask |= INITIAL_INTEGRATE;
    mask |= PRE_EXCHANGE;
    mask |= PRE_NEIGHBOR;
    mask |= PRE_FORCE;
    mask |= PRE_FINAL_INTEGRATE;
    mask |= FINAL_INTEGRATE;
    return mask;
}

/* ---------------------------------------------------------------------- */

double FixMultisphere::max_r_bound()
{
    return multisphere_.max_r_bound();
}

/* ---------------------------------------------------------------------- */

double FixMultisphere::extend_cut_ghost()
{

    return 2. * max_r_bound();
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::add_body_finalize()
{

    multisphere_.id_extend_body_extend(body_);
    multisphere_.generate_map();
    multisphere_.reset_forces(true);
    set_xv(LOOP_LOCAL);
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::init()
{
    // lots of error checks and warnings

    if (0 == atom->map_style)
        ms_error(FLERR, "requires an 'atom_modify map' command to allocate an atom map");

    if (!atom->rmass_flag || !atom->omega_flag)
        ms_error(FLERR, "need per-atom mass and omega");

    if (domain->dimension != 3)
        ms_error(FLERR, "works with 3D simulations only");

    if (modify->n_fixes_style("heat/gran") > 1)
        ms_error(FLERR, "only one fix heat/gran supported");
    fix_heat_ = static_cast<FixHeatGran *>(modify->find_fix_style("heat/gran", 0));

    if (fix_heat_ && atom->quaternion)
        ms_error(FLERR, "heat transfer not compatible with concave particles");

    if (domain->triclinic || dynamic_cast<DomainWedge *>(domain))
        ms_error(FLERR, "does not work with triclinic or wedge box");

    if (strstr(update->integrate_style, "respa"))
        ms_error(FLERR, "does not work with respa");
    // step_respa = ((Respa *) update->integrate)->step;

    if (force->newton)
        ms_error(FLERR, "requires newton 'off'");

    if (modify->n_fixes_style("gravity") > 1)
        ms_error(FLERR, "only one fix gravity supported");
    fix_gravity_ = static_cast<FixGravity *>(modify->find_fix_style("gravity", 0));

    // warn if more than one rigid fix
    if (modify->n_fixes_style("rigid") + modify->n_fixes_style("multisphere") > 1)
        error->warning(FLERR, "More than one fix rigid / fix multisphere");

    if (concave_)
    {
        int ntemp = modify->n_fixes_style("particletemplate");
        for (int itemp = 0; itemp < ntemp; itemp++)
        {
            Fix *tmp = modify->find_fix_style("particletemplate", itemp);
            if (strstr(tmp->style, "sphere"))
                ms_error(FLERR, "concave particles can not be combined with spheres, multispheres etc");
        }
    }

    fix_remove_.clear();

    // timestep info

    dtv = update->dt;
    dtf = 0.5 * update->dt * force->ftm2v;
    dtq = 0.5 * update->dt;
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::add_remove_callback(FixRemove *ptr)
{
    fix_remove_.push_back(ptr);

    if (atom->quaternion)
        ms_error(FLERR, "fix remove not compatible with concave particles");
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::setup(int vflag)
{
    // calc MS comm properties
    ntypes_ = modify->n_fixes_style("particletemplate/multisphere");

    ScalarContainer<int> *clumptypes = data().prop().getElementProperty<ScalarContainer<int>>("clumptype");
    int ntypes_existing = clumptypes->max();

    int nfmscfd = modify->n_fixes_style("couple/cfd/force/multisphere");

    if (ntypes_existing > ntypes_ && nfmscfd > 0)
        ms_error(FLERR, "for cfd coupling with multisphere drag force, you need to specify all "
                        "fix particletemplate/multisphere commands in case of restart that you had in the original set-up");

    if (Vclump_)
        delete[] Vclump_;
    Vclump_ = new double[ntypes_ + 1];

    for (int ifix = 0; ifix < ntypes_; ifix++)
    {
        FixTemplateMultisphere *ftm = static_cast<FixTemplateMultisphere *>(modify->find_fix_style("particletemplate/multisphere", ifix));
        int itype = ftm->type();
        Vclump_[itype] = ftm->volexpect();
    }

    int i, n;
    int nlocal = atom->nlocal;

    // virial setup before call to set_v

    if (vflag)
        v_setup(vflag);
    else
        evflag = 0;

    if (vflag_global)
        for (n = 0; n < 6; n++)
            virial[n] *= 2.0;
    if (vflag_atom)
    {
        for (i = 0; i < nlocal; i++)
            for (n = 0; n < 6; n++)
                vatom[i][n] *= 2.0;
    }

    if (fix_heat_ && !allow_heatsource_)
    {
        // check if heatsource is active for multisphere particles

        for (int i = 0; i < nlocal; i++)
        {
            // skip if atom not in rigid body
            if (body_[i] < 0)
                continue;

            int ibody = map(body_[i]);

            // skip if body not owned by this proc
            if (ibody < 0)
                continue;

            if (!domain->is_owned_or_first_ghost(i))
                continue;

            if (!MathExtraLiggghts::compDouble(fix_heat_->fix_heatSource->vector_atom[i], 0., 1e-6))
                ms_error(FLERR, "The multisphere heattransfer does not support heatsources");
        }
    }

    comm_correct_force(true);
    calc_force(true);
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::setup_pre_exchange()
{

    pre_exchange();
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::setup_pre_neighbor()
{

    pre_neighbor();
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::set_arrays(int i)
{

    body_[i] = -1;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixMultisphere::copy_arrays(int i, int j, int delflag)
{
    body_[j] = body_[i];
    displace_[j][0] = displace_[i][0];
    displace_[j][1] = displace_[i][1];
    displace_[j][2] = displace_[i][2];
}

/* ---------------------------------------------------------------------- */
/*
int FixMultisphere::getMask(int ibody)
{

  int    *mask     = atom->mask;
  int nloop = 0;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;

  nloop = nlocal+nghost;
  int counter = 0;
  int mask_curr = 0, mask_prev = 0;
  for (int i = 0; i < nloop; i++) {
    if (body_[i] < 0) continue;
    if(ibody == map(body_[i])) {
      mask_curr = mask[i];
      if(counter > 0) {
        if(mask_prev != mask_curr)
          error->one(FLERR,"Atoms in a multisphere particle have different group-IDs");
      }
      mask_prev = mask_curr;
      counter ++;
    }
  }
  return mask_curr;

}
*/

/* ---------------------------------------------------------------------- */

void FixMultisphere::initial_integrate(int vflag)
{
    int timestep = update->ntimestep;
    double **xcm = multisphere_.xcm_.begin();
    double **vcm = multisphere_.vcm_.begin();
    double **fcm = multisphere_.fcm_.begin();
    double **torquecm = multisphere_.torquecm_.begin();
    double **ex_space = multisphere_.ex_space_.begin();
    double **ey_space = multisphere_.ey_space_.begin();
    double **ez_space = multisphere_.ez_space_.begin();
    double **angmom = multisphere_.angmom_.begin();
    double **omega = multisphere_.omega_.begin();
    double **quat = multisphere_.quat_.begin();
    double **inertia = multisphere_.inertia_.begin();
    double *masstotal = multisphere_.masstotal_.begin();
    double *density = multisphere_.density_.begin();
    int *start_step = multisphere_.start_step_.begin();
    double **v_integrate = multisphere_.v_integrate_.begin();
    bool **fflag = multisphere_.fflag_.begin();
    bool **tflag = multisphere_.tflag_.begin();
    int nbody = multisphere_.n_body();

    if (strstr(style, "nointegration"))
        return;

    int n_stream = modify->n_fixes_style("insert/stream");
    bool has_stream = n_stream > 0;

    for (int ibody = 0; ibody < nbody; ibody++)
    {
        /*
        if(!(getMask(ibody) & groupbit) )
          continue;*/

        if (has_stream && timestep < start_step[ibody])
        {
            vectorCopy3D(v_integrate[ibody], vcm[ibody]);

            // update xcm by full step
            xcm[ibody][0] += dtv * vcm[ibody][0];
            xcm[ibody][1] += dtv * vcm[ibody][1];
            xcm[ibody][2] += dtv * vcm[ibody][2];

            continue;
        }

        // update vcm by 1/2 step

        const double addMassTerm = 1.0 + CAdd_ * fluidDensity_ / density[ibody];
        const double dtfm = dtf / (masstotal[ibody] * addMassTerm);

        if (fflag[ibody][0])
            vcm[ibody][0] += dtfm * fcm[ibody][0];
        if (fflag[ibody][1])
            vcm[ibody][1] += dtfm * fcm[ibody][1];
        if (fflag[ibody][2])
            vcm[ibody][2] += dtfm * fcm[ibody][2];

        // update xcm by full step

        xcm[ibody][0] += dtv * vcm[ibody][0];
        xcm[ibody][1] += dtv * vcm[ibody][1];
        xcm[ibody][2] += dtv * vcm[ibody][2];

        // update angular momentum by 1/2 step
        const double dtt = dtf / addMassTerm;
        if (tflag[ibody][0])
            angmom[ibody][0] += dtt * torquecm[ibody][0];
        if (tflag[ibody][1])
            angmom[ibody][1] += dtt * torquecm[ibody][1];
        if (tflag[ibody][2])
            angmom[ibody][2] += dtt * torquecm[ibody][2];

        // compute omega at 1/2 step from angmom at 1/2 step and current q
        // update quaternion a full step via Richardson iteration
        // returns new normalized quaternion, also updated omega at 1/2 step
        // update ex,ey,ez to reflect new quaternion

        MathExtra::angmom_to_omega(angmom[ibody], ex_space[ibody], ey_space[ibody],
                                   ez_space[ibody], inertia[ibody], omega[ibody]);
        MathExtra::richardson(quat[ibody], angmom[ibody], omega[ibody],
                              inertia[ibody], dtq);
        MathExtra::q_to_exyz(quat[ibody],
                             ex_space[ibody], ey_space[ibody], ez_space[ibody]);
    }

    // virial setup before call to set_xv

    if (vflag)
        v_setup(vflag);
    else
        evflag = 0;

    // set coords/orient and velocity/rotation of atoms in rigid bodies
    // from quarternion and omega

    set_xv();

    rev_comm_flag_ = MS_COMM_REV_X_V_OMEGA;
    reverse_comm();
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::setup_pre_force(int dummy)
{
    pre_force(dummy);
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::pre_force(int)
{
    // set force and torque to 0
    // do not reset external torques
    // other commands can use multisphere_.add_external_force()
    // in post_force

    multisphere_.reset_forces(false);
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::pre_final_integrate()
{
    comm_correct_force(false);
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::final_integrate()
{
    int timestep = update->ntimestep;
    // double **xcm = multisphere_.xcm_.begin();
    double **vcm = multisphere_.vcm_.begin();
    double **fcm = multisphere_.fcm_.begin();
    double **torquecm = multisphere_.torquecm_.begin();
    double **ex_space = multisphere_.ex_space_.begin();
    double **ey_space = multisphere_.ey_space_.begin();
    double **ez_space = multisphere_.ez_space_.begin();
    double **angmom = multisphere_.angmom_.begin();
    double **omega = multisphere_.omega_.begin();
    double **inertia = multisphere_.inertia_.begin();
    double *masstotal = multisphere_.masstotal_.begin();
    double *density = multisphere_.density_.begin();
    int *start_step = multisphere_.start_step_.begin();
    bool **fflag = multisphere_.fflag_.begin();
    bool **tflag = multisphere_.tflag_.begin();
    int nbody = multisphere_.n_body();

    // calculate forces and torques on body

    calc_force(false);

    if (strstr(style, "nointegration"))
        return;

    int n_stream = modify->n_fixes_style("insert/stream");
    bool has_stream = n_stream > 0;

    // resume integration
    for (int ibody = 0; ibody < nbody; ibody++)
    {
        /*
        if (!(getMask(ibody) & groupbit))
          continue; */

        if (has_stream && timestep < start_step[ibody])
            continue;

        // update vcm by 1/2 step
        const double addMassTerm = 1.0 + CAdd_ * fluidDensity_ / density[ibody];
        const double dtfm = dtf / (masstotal[ibody] * addMassTerm);
        if (fflag[ibody][0])
            vcm[ibody][0] += dtfm * fcm[ibody][0];
        if (fflag[ibody][1])
            vcm[ibody][1] += dtfm * fcm[ibody][1];
        if (fflag[ibody][2])
            vcm[ibody][2] += dtfm * fcm[ibody][2];

        // update angular momentum by 1/2 step
        const double dtt = dtf / addMassTerm;
        if (tflag[ibody][0])
            angmom[ibody][0] += dtt * torquecm[ibody][0];
        if (tflag[ibody][1])
            angmom[ibody][1] += dtt * torquecm[ibody][1];
        if (tflag[ibody][2])
            angmom[ibody][2] += dtt * torquecm[ibody][2];

        MathExtra::angmom_to_omega(angmom[ibody], ex_space[ibody], ey_space[ibody],
                                   ez_space[ibody], inertia[ibody], omega[ibody]);
    }

    set_v();

    rev_comm_flag_ = MS_COMM_REV_V_OMEGA;
    reverse_comm();

    fw_comm_flag_ = MS_COMM_FW_V_OMEGA;
    forward_comm();
}

/* ----------------------------------------------------------------------
   call to set_v, plus according fwd communication
------------------------------------------------------------------------- */

void FixMultisphere::set_v_communicate()
{

    set_v();

    rev_comm_flag_ = MS_COMM_REV_V_OMEGA;
    reverse_comm();

    fw_comm_flag_ = MS_COMM_FW_V_OMEGA;
    forward_comm();
}

/* ----------------------------------------------------------------------
   communicate and correct forces for multisphere bodies
------------------------------------------------------------------------- */

void FixMultisphere::comm_correct_force(bool setupflag)
{
    // communication and correction before real integration

    fw_comm_flag_ = MS_COMM_FW_F_TORQUE;
    forward_comm();

    if (setupflag)
        fix_volumeweight_ms_->do_forward_comm();

    // correct forces if necessary
    if (do_modify_body_forces_torques_)
        modify_body_forces_torques();
}

/* ----------------------------------------------------------------------
   set space-frame coords and velocity of each atom in each rigid body
   set orientation and rotation of extended particles
   x = Q displace + Xcm, mapped back to periodic box
   v = Vcm + (W cross (x - Xcm))
------------------------------------------------------------------------- */

void FixMultisphere::calc_force(bool setupflag)
{
    int ibody;
    tagint *image = atom->image;
    double **x = atom->x;
    double **f_atom = atom->f;
    double **torque_atom = atom->torque;
    double f_one[3], torque_one[3];
    int nlocal = atom->nlocal;
    int nghost = atom->nghost;

    double **xcm = multisphere_.xcm_.begin();
    double *masstotal = multisphere_.masstotal_.begin();
    double **fcm = multisphere_.fcm_.begin();
    double **dragforce_cm = multisphere_.dragforce_cm_.begin();
    double **hdtorque_cm = multisphere_.hdtorque_cm_.begin();
    double **torquecm = multisphere_.torquecm_.begin();
    double *temp = multisphere_.temp_.begin();
    double *temp_old = multisphere_.temp_old_.begin();
    int nbody = multisphere_.n_body();

    double unwrap[3], dx, dy, dz;

    // calculate forces and torques of bodies
    for (int i = 0; i < nlocal + nghost; i++)
    {   
        if (domain->lebc)
        {
            force_length[i][0] = 0.0;
            force_length[i][1] = 0.0;
            force_length[i][2] = 0.0;
            force_length[i][3] = 0.0;
            force_length[i][4] = 0.0;
            force_length[i][5] = 0.0;
        }

        if (body_[i] < 0)
            continue;

        ibody = map(body_[i]);

        if (ibody < 0)
            continue;

        if (!domain->is_owned_or_first_ghost(i))
            continue;

        vectorCopy3D(f_atom[i], f_one);
        vectorCopy3D(torque_atom[i], torque_one);

        fcm[ibody][0] += f_one[0];
        fcm[ibody][1] += f_one[1];
        fcm[ibody][2] += f_one[2];

        domain->unmap(x[i], image[i], unwrap); // TODO?!
        dx = unwrap[0] - xcm[ibody][0];
        dy = unwrap[1] - xcm[ibody][1];
        dz = unwrap[2] - xcm[ibody][2];

        if (i >= nlocal)
            domain->minimum_image(dx, dy, dz); // TODO?!

        torquecm[ibody][0] += dy * f_one[2] - dz * f_one[1] + torque_one[0];
        torquecm[ibody][1] += dz * f_one[0] - dx * f_one[2] + torque_one[1];
        torquecm[ibody][2] += dx * f_one[1] - dy * f_one[0] + torque_one[2];

        if (domain->lebc)
        {
            force_length[i][0] =   -dx * f_one[0];
            force_length[i][1] =   -dy * f_one[1];
            force_length[i][2] =   -dz * f_one[2];
            force_length[i][3] =   -dx * f_one[1];
            force_length[i][4] =   -dx * f_one[2];
            force_length[i][5] =   -dy * f_one[2];

        }
    }

    // heat transfer

    if (fix_heat_)
    {
        // communicate temperature to ghosts
        fw_comm_flag_ = MS_COMM_FW_TEMP;
        forward_comm();

        // save old temp
        for (ibody = 0; ibody < nbody; ibody++)
            temp_old[ibody] = temp[ibody];

        if (setupflag)
            vectorZeroizeN(temp, nbody);

        // caclulate temperature from single particles
        for (int i = 0; i < nlocal + nghost; i++)
        {
            // skip if atom not in rigid body
            if (body_[i] < 0)
                continue;

            ibody = map(body_[i]);

            // skip if body not owned by this proc
            if (ibody < 0)
                continue;

            if (!domain->is_owned_or_first_ghost(i))
                continue;

            // if (screen) fprintf(screen, "Update temperature for particle i = %d\n",i);
            // if (screen) fprintf(screen, "Heat flux of particle %d is %g\n",i,fix_heat_->heatFlux[i]);

            if (!setupflag)
                temp[ibody] += fix_heat_->fix_temp->vector_atom[i] - temp_old[ibody]; // fix_heat_->heatFlux[i]*update->dt/(masstotal[ibody]);
            else
                temp[ibody] += fix_volumeweight_ms_->vector_atom[i] * fix_heat_->fix_temp->vector_atom[i];
        }

        // set temperature of single particles
        for (int i = 0; i < nlocal + nghost; i++)
        {
            // skip if atom not in rigid body
            if (body_[i] < 0)
                continue;

            ibody = map(body_[i]);

            // skip if body not owned by this proc
            if (ibody < 0)
                continue;

            if (!domain->is_owned_or_first_ghost(i))
                continue;

            fix_heat_->fix_temp->vector_atom[i] = temp[ibody];
        }

        rev_comm_flag_ = MS_COMM_REV_TEMP;
        reverse_comm();
    }

    // add external forces on bodies, such as gravity, dragforce

    if (fix_gravity_)
    {
        double grav[3];
        fix_gravity_->get_gravity(grav);
        for (ibody = 0; ibody < nbody; ibody++)
        {
            fcm[ibody][0] += masstotal[ibody] * grav[0];
            fcm[ibody][1] += masstotal[ibody] * grav[1];
            fcm[ibody][2] += masstotal[ibody] * grav[2];
        }
    }
    if (add_dragforce_)
    {
        for (ibody = 0; ibody < nbody; ibody++)
        {

            vectorAdd3D(fcm[ibody], dragforce_cm[ibody], fcm[ibody]);
            vectorAdd3D(torquecm[ibody], hdtorque_cm[ibody], torquecm[ibody]);
        }
    }
}

/* ----------------------------------------------------------------------
   set space-frame coords and velocity of each atom in each rigid body
   set orientation and rotation of extended particles
   x = Q displace + Xcm, mapped back to periodic box
   v = Vcm + (W cross (x - Xcm))
------------------------------------------------------------------------- */

void FixMultisphere::set_xv()
{
    set_xv(LOOP_ALL);
}

/*-----------------------------------------------------------------------*/

void FixMultisphere::set_xv(int ghostflag)
{
    int ibody;
    int xbox, ybox, zbox;
    double x0 = 0.0, x1 = 0.0, x2 = 0.0, v0 = 0.0, v1 = 0.0, v2 = 0.0, massone;
    double vr[6];

    tagint *image = atom->image;
    double **x = atom->x;
    double **v = atom->v;
    double **f = atom->f;
    double **omega_one = atom->omega;
    double *rmass = atom->rmass;
    double *mass = atom->mass;
    int *type = atom->type;
    int nlocal = atom->nlocal;
    int nghost = atom->nghost;

    double **xcm = multisphere_.xcm_.begin();
    double **vcm = multisphere_.vcm_.begin();
    double **ex_space = multisphere_.ex_space_.begin();
    double **ey_space = multisphere_.ey_space_.begin();
    double **ez_space = multisphere_.ez_space_.begin();
    double **omega = multisphere_.omega_.begin();
    double **quat = multisphere_.quat_.begin();

    int nloop = 0;

    if (ghostflag == LOOP_ALL)
        nloop = nlocal + nghost;
    else if (ghostflag == LOOP_LOCAL)
        nloop = nlocal;
    else
        error->all(FLERR, "Illegal call to FixMultisphere::set_v");

    double xprd = domain->xprd;
    double yprd = domain->yprd;
    double zprd = domain->zprd;

    // set x and v of each atom

    for (int i = 0; i < nloop; i++)
    {

        if (body_[i] < 0)
            continue;
        ibody = map(body_[i]);

        if (ibody < 0)
            continue;

        xbox = (image[i] & IMGMASK) - IMGMAX;
        ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (image[i] >> IMG2BITS) - IMGMAX;

        // save old positions and velocities for virial

        if (evflag)
        {
            x0 = x[i][0] + xbox * xprd + ybox * domain->lebc_displacement;
            x1 = x[i][1] + ybox * yprd;
            x2 = x[i][2] + zbox * zprd;
            v0 = v[i][0] + ybox * domain->yprd * domain->ssr;
            v1 = v[i][1];
            v2 = v[i][2];
        }

        // x = displacement from center-of-mass, based on body orientation
        // v = vcm + omega around center-of-mass

        MathExtra::matvec(ex_space[ibody], ey_space[ibody], ez_space[ibody], displace_[i], x[i]);
        if (body_[i] == 3)
        {
            // if (i >= nlocal)
            // {
            //     if (ybox != 0 || xbox != 0)
            //         std::cout << xbox << " " << ybox << " ghost\n";
            // }
            if (i < nlocal)
            {
                if (ybox != 0 || xbox != 0)
                {
                    
                    int values[4];
                    data().remap_values(values, ibody);
                    if (values[0] != 0 ||values[1] != 0 ||values[2] != 0 ||values[3] != 0)
                    {
                        std::cout << xbox << " " << ybox << " local " << comm->me << "\n" << std::flush;
                        std::cout << "remap " << values[0] << values[1] << values[2] << values[3] << "\n" << std::flush;
                    }
                        

                }
                
                
            }
            
        }
        

        v[i][0] = omega[ibody][1] * x[i][2] - omega[ibody][2] * x[i][1] + vcm[ibody][0] - ybox * domain->yprd * domain->ssr; // TODO?!
        v[i][1] = omega[ibody][2] * x[i][0] - omega[ibody][0] * x[i][2] + vcm[ibody][1];
        v[i][2] = omega[ibody][0] * x[i][1] - omega[ibody][1] * x[i][0] + vcm[ibody][2];

        // add center of mass to displacement
        // map back into periodic box via xbox,ybox,zbox
        // for triclinic, would have to add in box tilt factors as well

        x[i][0] += xcm[ibody][0] - xbox * xprd - ybox * domain->lebc_displacement;
        x[i][1] += xcm[ibody][1] - ybox * yprd;
        x[i][2] += xcm[ibody][2] - zbox * zprd;

        omega_one[i][0] = omega[ibody][0];
        omega_one[i][1] = omega[ibody][1];
        omega_one[i][2] = omega[ibody][2];

        // set quat as well if applicable
        if (atom->quaternion)
            vectorCopy4D(quat[ibody], atom->quaternion[i]);

        // virial = unwrapped coords dotted into body constraint force
        // body constraint force = implied force due to v change minus f external
        // assume f does not include forces internal to body
        // 1/2 factor b/c final_integrate contributes other half
        // assume per-atom contribution is due to constraint force on that atom

        if (evflag && i < nlocal)
        {
            if (rmass)
                massone = rmass[i];
            else
                massone = mass[type[i]];
            const double fc0 = massone * (v[i][0] - v0) / dtf - f[i][0];
            const double fc1 = massone * (v[i][1] - v1) / dtf - f[i][1];
            const double fc2 = massone * (v[i][2] - v2) / dtf - f[i][2];

            vr[0] = 0.5 * x0 * fc0;
            vr[1] = 0.5 * x1 * fc1;
            vr[2] = 0.5 * x2 * fc2;
            vr[3] = 0.5 * x0 * fc1;
            vr[4] = 0.5 * x0 * fc2;
            vr[5] = 0.5 * x1 * fc2;

            v_tally(1, &i, 1.0, vr);
        }
    }
}

/* ----------------------------------------------------------------------
   set space-frame velocity of each atom in a rigid body
   set omega and angmom of extended particles
   v = Vcm + (W cross (x - Xcm))
------------------------------------------------------------------------- */

void FixMultisphere::set_v()
{
    set_v(LOOP_ALL);
}

/*-----------------------------------------------------------------------*/

void FixMultisphere::set_v(int ghostflag)
{
    int ibody;
    int xbox, ybox, zbox;
    double x0 = 0.0, x1 = 0.0, x2 = 0.0, v0 = 0.0, v1 = 0.0, v2 = 0.0, massone;
    double delta[3], vr[6];

    double **x = atom->x;
    double **v = atom->v;
    double **f = atom->f;
    double *rmass = atom->rmass;
    double *mass = atom->mass;
    double **omega_one = atom->omega;
    int *type = atom->type;
    tagint *image = atom->image;
    int nlocal = atom->nlocal;
    int nghost = atom->nghost;
    double **vcm = multisphere_.vcm_.begin();
    double **omega = multisphere_.omega_.begin();
    double **ex_space = multisphere_.ex_space_.begin();
    double **ey_space = multisphere_.ey_space_.begin();
    double **ez_space = multisphere_.ez_space_.begin();

    int nloop = 0;

    if (ghostflag == LOOP_ALL)
        nloop = nlocal + nghost;
    else if (ghostflag == LOOP_LOCAL)
        nloop = nlocal;
    else
        error->all(FLERR, "Illegal call to FixMultisphere::set_v");

    double xprd = domain->xprd;
    double yprd = domain->yprd;
    double zprd = domain->zprd;

    // set v of each atom

    for (int i = 0; i < nloop; i++)
    {
        if (body_[i] < 0)
            continue;
        ibody = map(body_[i]);
        if (ibody < 0)
            continue;

        MathExtra::matvec(ex_space[ibody], ey_space[ibody], ez_space[ibody], displace_[i], delta);

        xbox = (image[i] & IMGMASK) - IMGMAX;
        ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (image[i] >> IMG2BITS) - IMGMAX;

        // save old velocities for virial

        if (evflag)
        {
            v0 = v[i][0]; // + ybox * domain->yprd * domain->ssr;
            v1 = v[i][1];
            v2 = v[i][2];
        }

        // if (i >= nlocal)
        // {
        //     if (ybox != 0 && xbox != 0)
        //         std::cout << xbox << " " << ybox << "\n";
        // }

        v[i][0] = omega[ibody][1] * delta[2] - omega[ibody][2] * delta[1] + vcm[ibody][0] - ybox * domain->yprd * domain->ssr; // TODO?!
        v[i][1] = omega[ibody][2] * delta[0] - omega[ibody][0] * delta[2] + vcm[ibody][1];
        v[i][2] = omega[ibody][0] * delta[1] - omega[ibody][1] * delta[0] + vcm[ibody][2];

        omega_one[i][0] = omega[ibody][0];
        omega_one[i][1] = omega[ibody][1];
        omega_one[i][2] = omega[ibody][2];

        // virial = unwrapped coords dotted into body constraint force
        // body constraint force = implied force due to v change minus f external
        // assume f does not include forces internal to body
        // 1/2 factor b/c initial_integrate contributes other half
        // assume per-atom contribution is due to constraint force on that atom

        if (evflag && i < nlocal)
        {
            if (rmass)
                massone = rmass[i];
            else
                massone = mass[type[i]];
            const double fc0 = massone * (v[i][0] - v0) / dtf - f[i][0];
            const double fc1 = massone * (v[i][1] - v1) / dtf - f[i][1];
            const double fc2 = massone * (v[i][2] - v2) / dtf - f[i][2];

            x0 = x[i][0] + xbox * xprd + ybox * domain->lebc_displacement;
            x1 = x[i][1] + ybox * yprd;
            x2 = x[i][2] + zbox * zprd;

            vr[0] = 0.5 * x0 * fc0;
            vr[1] = 0.5 * x1 * fc1;
            vr[2] = 0.5 * x2 * fc2;
            vr[3] = 0.5 * x0 * fc1;
            vr[4] = 0.5 * x0 * fc2;
            vr[5] = 0.5 * x1 * fc2;

            v_tally(1, &i, 1.0, vr);
        }
    }
}

/* ----------------------------------------------------------------------
   delete atoms belonging to deleted bodies
------------------------------------------------------------------------- */

void FixMultisphere::pre_exchange()
{
    AtomVec *avec = atom->avec;

    // reset last trigger for re-neigh
    next_reneighbor = -1;

    double *delflag = fix_delflag_->vector_atom;

    int i = 0;

    while (i < atom->nlocal)
    {

        if (MathExtraLiggghts::compDouble(delflag[i], 1.))
        {

            avec->copy(atom->nlocal - 1, i, 1);
            atom->nlocal--;
        }
        else
            i++;
    }

    
    if (force_volume_fraction)
    {
        pti_list = (ParticleToInsertMultisphere **)memory->smalloc(data().deleted_body_types.size() * sizeof(ParticleToInsertMultisphere *), "pti_list");
        int count = 0;
        
        for (int body_type : data().deleted_body_types)
        {
            std::cout << "fix_ms forced to insert " << body_type << " body types\n";
            FixTemplateMultisphere *fix_ms_template = static_cast<FixTemplateMultisphere *>(modify->find_fix_style("particletemplate/multisphere", body_type-1));

            if (fix_ms_template)
            {   
                if (fix_ms_template->type() != body_type)
                {
                    std::cout << "fix_ms template wrong body types\n";
                }
                
                pti_list[count] = fix_ms_template->get_particle_to_insert();
               
                pti_list[count]->set_x_v_omega_center_processor(count, data().deleted_body_types.size());
              
                
                count++;
            } else {
                std::cout << "not found: fix_ms template\n";
            }
        }
        int n = 0;
        for (int i = 0; i < count; i++)
        {   
           
            n += pti_list[i]->insert();
        }
      
        atom->tag_extend();
        bigint all_atoms;
        bigint added_atoms = n;
        MPI_Allreduce(&added_atoms, &all_atoms, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
       
        atom->natoms += all_atoms;// why just add here?

        atom->nghost = 0;
        atom->map_init();
        atom->map_set();
       

        add_body_finalize();
        for (int i = 0; i < data().deleted_body_types.size(); i++)
        {
            delete pti_list[i];
        }

        data().deleted_body_types.clear();
        memory->sfree(pti_list);
     

    }
   
}

/* ----------------------------------------------------------------------
   communicate body

   remap xcm of each rigid body back into periodic simulation box
   done during pre_neighbor so will be after call to pbc()
     and after fix_deform::pre_exchange() may have flipped box
   use domain->remap() in case xcm is far away from box
     due to 1st definition of rigid body or due to box flip
   if don't do this, then atoms of a body which drifts far away
     from a triclinic box will be remapped back into box
     with huge displacements when the box tilt changes via set_x()

   exchange bodies with stencil procs

   check for lost bodies and remove them
     also mark atoms belonging to lost bodies for deletion

   communicate displace, image
------------------------------------------------------------------------- */

void FixMultisphere::pre_neighbor()
{

    int nall = atom->nlocal + atom->nghost;
    double *corner_ghost = fix_corner_ghost_->vector_atom;
    vectorZeroizeN(corner_ghost, nall);



    fw_comm_flag_ = MS_COMM_FW_BODY;
    forward_comm();



    for (size_t irem = 0; irem < fix_remove_.size(); irem++)
        (fix_remove_[irem])->delete_bodies();



    fw_comm_flag_ = MS_COMM_FW_IMAGE_DISPLACE;
    forward_comm();



    multisphere_.remap_bodies(body_);


    rev_comm_flag_ = MS_COMM_REV_IMAGE;
    reverse_comm();


    multisphere_.exchange();

    
    multisphere_.calc_nbody_all();

    if (lebc_highest_body_id != multisphere_.n_body_all() && domain->lebc)
    {
        std::cout << "bodies changed from: " << lebc_highest_body_id << " to: " << multisphere_.n_body_all() << "\n"<<std::flush;
        MPI_Abort(0,1);
    }

    multisphere_.generate_map();

    // set deletion flag
    // if any deleted atoms, do re-neigh in 100 steps at latest to remove
    // remainder particles
    double *delflag = fix_delflag_->vector_atom;
    double *existflag = fix_existflag_->vector_atom;
    vectorZeroizeN(delflag, atom->nlocal + atom->nghost);
    vectorZeroizeN(existflag, atom->nlocal + atom->nghost);

    bigint tmp = atom->nlocal;
    MPI_Allreduce(&tmp, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);


    if(check_bodies_atoms("FixMultisphere::pre_neighbor()"))
    {
        std::cout << atom->natoms << " total atoms\n" <<std::flush;
    }



    if (multisphere_.check_lost_atoms(body_, delflag, existflag, fix_volumeweight_ms_->vector_atom))
        next_reneighbor = update->ntimestep + 5;

    fix_delflag_->do_reverse_comm();
    fix_existflag_->do_reverse_comm();

    fw_comm_flag_ = MS_COMM_FW_IMAGE_DISPLACE;
    forward_comm();

    // merge delflag and existflag

    int nlocal = atom->nlocal;
    delflag = fix_delflag_->vector_atom;
    existflag = fix_existflag_->vector_atom;
    int forceNeighbour = 0; // use int instead of bool for MPI
    for (int i = 0; i < nlocal; i++)
    {
        // fprintf(screen,"On proc %d: For ilocal = %d is existflag = %g and delflag = %g \n",comm->me,i,existflag[i],delflag[i]);
        delflag[i] = (MathExtraLiggghts::compDouble(existflag[i], 0.)) ? 1. : delflag[i];
        if (MathExtraLiggghts::compDouble(delflag[i], 1.0))
            forceNeighbour = 1;
    }
    MPI_Max_Scalar(forceNeighbour, world);
    if (forceNeighbour)
        next_reneighbor = update->ntimestep + 5;
    

    multisphere_.calc_nbody_all();

    if (lebc_highest_body_id != multisphere_.n_body_all() && domain->lebc)
    {
        std::cout << "at end bodies changed from: " << lebc_highest_body_id << " to: " << multisphere_.n_body_all() << "\n" <<std::flush;
        MPI_Abort(0,1);
    }
}

/* ----------------------------------------------------------------------
   count # of degrees-of-freedom removed by fix_rigid for atoms in igroup
------------------------------------------------------------------------- */

int FixMultisphere::dof(int igroup)
{
    int n = 0;

    if (0 == comm->me)
        error->warning(FLERR, "Energy calculated for multisphere particles is currently not correct");
    // error->all(FLERR,"For multisphere particles, please use compute ke and instead of the thermo keyword ke");

    return n;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixMultisphere::memory_usage()
{
    int nmax = atom->nmax;
    double bytes = nmax * sizeof(int);
    bytes += nmax * 3 * sizeof(double);
    bytes += maxvatom * 6 * sizeof(double);

    // add Multisphere memory usage

    return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixMultisphere::grow_arrays(int nmax)
{

    body_ = memory->grow(body_, nmax, "rigid:body_");
    memory->grow(displace_, nmax, 3, "rigid:displace");
    atom->molecule = body_;

    memory->grow(force_length, nmax, 6, "rigid:force_length");

}

/* ----------------------------------------------------------------------
   extract values
------------------------------------------------------------------------- */
/*
void * FixMultisphere::extract(char *name, int &len1, int &len2)
{
    return multisphere_.extract(name,len1,len2);
}*/

/* ----------------------------------------------------------------------
   return attributes of a rigid body
   12 values per body
   xcm = 1,2,3; vcm = 4,5,6; fcm = 7,8,9; torque = 10,11,12
------------------------------------------------------------------------- */

double **FixMultisphere::get_dump_ref(int &nb, int &nprop, char *prop)
{
    error->one(FLERR, "TODO");
    return NULL;
}

/* ----------------------------------------------------------------------
   size and maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixMultisphere::size_restart(int nlocal)
{
    return 4 + 1;
}

int FixMultisphere::maxsize_restart()
{
    return 4 + 1;
}

/* ----------------------------------------------------------------------
   restart global
------------------------------------------------------------------------- */

void FixMultisphere::write_restart(FILE *fp)
{
    multisphere_.writeRestart(fp);
}

void FixMultisphere::restart(char *buf)
{
    double *list = (double *)buf;

    bool have_massflow_mesh = modify->have_restart_data_style("massflow/mesh");
    if (have_massflow_mesh)
    {
        int nmassflow = modify->n_restart_data_global_style("massflow/mesh");

        char property_name[200];
        for (int imf = 0; imf < nmassflow; imf++)
        {
            char *id_this = modify->id_restart_data_global_style("massflow/mesh", imf);
            sprintf(property_name, "counter_ms_%s", id_this);
            multisphere_.prop().addElementProperty<ScalarContainer<int>>(static_cast<const char *>(property_name), "comm_exchange_borders", "frame_invariant", "restart_yes");
        }
    }

    multisphere_.restart(list);
}

/* ---------------------------------------------------------------------- */

int FixMultisphere::modify_param(int narg, char **arg)
{
    if (strcmp(arg[0], "fflag") == 0)
    {
        bool fflag[3] = {true, true, true};
        if (narg < 4)
            ms_error(FLERR, "not enough arguments for 'fflag'");
        if (strcmp(arg[1], "on") == 0)
            fflag[0] = true;
        else if (strcmp(arg[1], "off") == 0)
            fflag[0] = false;
        else
            ms_error(FLERR, "expecting 'on or 'off' after 'fflag'");
        if (strcmp(arg[2], "on") == 0)
            fflag[1] = true;
        else if (strcmp(arg[2], "off") == 0)
            fflag[1] = false;
        else
            ms_error(FLERR, "expecting 'on or 'off' after 'fflag'");
        if (strcmp(arg[3], "on") == 0)
            fflag[2] = true;
        else if (strcmp(arg[3], "off") == 0)
            fflag[2] = false;
        else
            ms_error(FLERR, "expecting 'on or 'off' after 'fflag'");

        // MODIFY FFLAGS FOR ALL MS
        int nbody = multisphere_.n_body();
        for (int ibody = 0; ibody < nbody; ibody++)
        {
            multisphere_.set_fflag(ibody, fflag);
        }
        return 4;
    }
    else if (strcmp(arg[0], "tflag") == 0)
    {
        bool tflag[3] = {true, true, true};
        if (narg < 4)
            ms_error(FLERR, "not enough arguments for 'tflag'");
        if (strcmp(arg[1], "on") == 0)
            tflag[0] = true;
        else if (strcmp(arg[1], "off") == 0)
            tflag[0] = false;
        else
            ms_error(FLERR, "expecting 'on or 'off' after 'tflag'");
        if (strcmp(arg[2], "on") == 0)
            tflag[1] = true;
        else if (strcmp(arg[2], "off") == 0)
            tflag[1] = false;
        else
            ms_error(FLERR, "expecting 'on or 'off' after 'tflag'");
        if (strcmp(arg[3], "on") == 0)
            tflag[2] = true;
        else if (strcmp(arg[3], "off") == 0)
            tflag[2] = false;
        else
            ms_error(FLERR, "expecting 'on or 'off' after 'tflag'");

        // MODIFY TFLAGS FOR ALL MS
        int nbody = multisphere_.n_body();
        for (int ibody = 0; ibody < nbody; ibody++)
        {
            multisphere_.set_tflag(ibody, tflag);
        }
        return 4;
    }
    return 0;
}



bool FixMultisphere::check_bodies_atoms(std::string fff){

//   return;
  if (!domain->lebc) return false;
  
  
  std::vector<int> error_bodies = multisphere_.check_lost_atoms_no_delete(body_);

  bool found_error = false;
  for (int lob: error_bodies) 
  {
    if (lob >= 0)
    {
        // std::cout << lob << "\n";
        if (domain->lebc)
        {   
            found_error = true;
            std::cout <<std::to_string(update->ntimestep)<<" "<<fff<<", lost_body_tag="<<lob<<"\n"<<std::flush;

            for (int i = 0; i < atom->nlocal + atom->nghost; i++) {
                int body_tag_from_atom = body_[i];

                char real_or_ghost = 'r';
                if (i >= atom->nlocal)
                {
                    real_or_ghost = 'g';
                }
                

                if (body_tag_from_atom==lob) {
                std::cout <<std::to_string(update->ntimestep)
                            <<" "<<fff<< ": "<<real_or_ghost<< " atom_tag="<<atom->tag[i] << " atom_map="<<atom->map(atom->tag[i])
                            <<", x/L="<<atom->x[i][0]/domain->prd[0]
                            <<", y/L="<<atom->x[i][1]/domain->prd[1]
                            <<", z/L="<<atom->x[i][2]/domain->prd[2]
                            <<", bt="<<belongs_to(i)
                            <<", me="<< comm->me
                            <<"\n"<<std::flush;
                }
            }
        }
    }
  }
  return found_error;
//   int print = 1;

//   tagint *image0 = atom->image;
//   std::cout<<std::scientific;
//   int ex = -1;
//   double ssr = domain->ssr;
//   double dv = ssr * domain->prd[1];
//   double maxv = 0.0;
//   for (int i = 0; i < n_body(); i++) {
//     if (tag(i) > lebc_highest_body_id) {
//       std::vector<int> iatoms = atoms_of_body(tag(i), atom->nlocal + atom->nghost);
//       std::vector<int> imgs(3); 
//       int j = 0;
//       double v0 = 0.0,v1=0.0;
//       for (int iatom : iatoms)  {
//         int xbox = (image0[iatom] & IMGMASK) - IMGMAX;
//         int ybox = (image0[iatom] >> IMGBITS & IMGMASK) - IMGMAX;
//         int zbox = (image0[iatom] >> IMG2BITS) - IMGMAX;
//         if (j==0) {
//           v0 = std::sqrt(atom->v[iatom][0]*atom->v[iatom][0]+
//                          atom->v[iatom][1]*atom->v[iatom][1]+
//                          atom->v[iatom][2]*atom->v[iatom][2]);
//           imgs[0] = xbox;
//           imgs[1] = ybox;
//           imgs[2] = zbox;
//           j++;
//         } else {
//           v1 = std::sqrt(atom->v[iatom][0]*atom->v[iatom][0]+
//                          atom->v[iatom][1]*atom->v[iatom][1]+
//                          atom->v[iatom][2]*atom->v[iatom][2]);
//           // check x box images
//           if (imgs[0]!=xbox && (v1+v0)/2.0>dv) {
//             maxv = (v1+v0)/2.0;
//             // found one body worth looking into
//             ex = tag(i);
//             std::cout <<std::to_string(update->ntimestep)<<" "<<fff<< ": body_id="<<ex<<", maxv="<<maxv<<", dv="<<dv<<"\n"<<std::flush;
//             break;
//           }
//         }
//       }
//     }
//     if (ex!=-1) break;
//   }
//   return;
  // // at this point, if one of the cores has a body to investigate:
  // double gmaxv;
  // MPI_Allreduce(&maxv,&gmaxv,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  // // check who has this max
  // int sroot = -1;
  // if (maxv==gmaxv) {
  //   sroot = comm->me;
  // }
  // MPI_Barrier(MPI_COMM_WORLD);
  // MPI_Bcast(&ex,1,MPI_DOUBLE,sroot,MPI_COMM_WORLD);

  // if (comm->me==0 && domain->switchLEbc && ex>0) std::cout <<std::to_string(update->ntimestep)<<" "<<fff<< ": body_id="<<ex<<"\n"<<std::flush;

  // if (ex<0) return;
  // double **v = atom->v;
  // double tmpx[3],tmpv[3];
  // tagint *image = atom->image;
  // double vmax_b[3] = {0};
  // double vmax_a[3] = {0};
  // double vmax_a2[3] = {0};
  // double vmax_avel = 0;
  // int nlocal = atom->nlocal;// + atom->nghost;
  // int aid = -1;
  // int aidp = -1;

  // for (int i = 0; i < nlocal; i++) {    
  //   double vel = sqrt(v[i][0]*v[i][0] + 
  // 		      v[i][1]*v[i][1] +
  // 		      v[i][2]*v[i][2]);
  //   if (vel>vmax_avel) vmax_avel = vel;    
  // }

  // double vmm = 0;
  // MPI_Allreduce(&vmax_avel, &vmm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);  
  
  // int comm_found = -1;
  // for (int i = 0; i < nlocal; i++) {    
  //   double vel = sqrt(v[i][0]*v[i][0] + 
  // 		      v[i][1]*v[i][1] +
  // 		      v[i][2]*v[i][2]);
  //   if (abs(vmm-vel)>1e-16){
  //     vmax_avel = vel; 
  //     vmax_a2[0] = v[i][0];
  //     vmax_a2[1] = v[i][1];
  //     vmax_a2[2] = v[i][2];
  //     aid=i;
  //     comm_found = comm->me;
  //   }
  // }
  // // return;
  // if (comm_found == comm->me) {
  //   aidp = aid;
  //   int body_index = multisphere_.map(belongs_to(aidp));

  //   const int binproc = (body_index>=0);
  //   if (belongs_to(aidp) <= 0 || !binproc) return;
  //   // if (fix_ms->tag(i) > fix_ms->lebc_highest_body_id || fix_ms->tag(i) <= 0  // || fix_ms->tag(i) == 19
  //   // 	) continue;

  //   get_vcm(tmpv, body_index);

  //   double tv[3]; 
  //   get_vcm(tv, body_index);
  //   tmpv[0] = abs((tmpv[0] - v[aidp][0])/ (abs(tmpv[0])>1e-16?tmpv[0]:1.0));
  //   tmpv[1] = abs((tmpv[1] - v[aidp][1])/ (abs(tmpv[1])>1e-16?tmpv[1]:1.0));
  //   tmpv[2] = abs((tmpv[2] - v[aidp][2])/ (abs(tmpv[2])>1e-16?tmpv[2]:1.0));
    
  //   if ((tmpv[0]>2.0 ||
  // 	 tmpv[1]>2.0 || 
  // 	 tmpv[2]>2.0 ) &&
  // 	(abs(v[aidp][0])>1 || abs(v[aidp][1])>1 || abs(v[aidp][2])>1)
  // 	)
  //     std::cout <<std::to_string(update->ntimestep)
  // 		<<" "<<fff<<", tag="<<tag(body_index)
  // 		<<", body=["<<tv[0]<<","<<tv[1]<<","<<tv[2]<<"]"
  // 		<<", diff=["<<tmpv[0]<<","<<tmpv[1]<<","<<tmpv[2]<<"]"
  // 		<<", atom=["<<v[aidp][0]<<","<<v[aidp][1]<<","<<v[aidp][2]<<"]\n";
    
    
  // }

  // // double vmax_b2 = 0;  
  // // if (aidp==aid) {
  // //   int body_index = fix_ms->data().map(fix_ms->belongs_to(aidp));
  // //   get_vcm(tmpv, body_index);
  // //   vmax_b2[p] = tmpv[p];
  // // }

  // // MPI_Reduce(&vmax_b2[0], &vmax_b[0], 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);  

  // // if ()
  // // int check = 0;
  // // if (abs(vmax_b[0])>1e-16) check = std::max(int(vmax_a[0]/vmax_b[0]>2.0  || vmax_a[0]/vmax_b[0]<0.5),check);
  // // if (abs(vmax_b[1])>1e-16) check = std::max(int(vmax_a[1]/vmax_b[1]>2.0  || vmax_a[1]/vmax_b[1]<0.5),check);
  // // if (abs(vmax_b[2])>1e-16) check = std::max(int(vmax_a[2]/vmax_b[2]>2.0  || vmax_a[2]/vmax_b[2]<0.5),check);

  // // int c2 = 0;
  // // MPI_Reduce(&check, &c2, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

  // // if (comm->me==0 && c2) std::cout <<std::to_string(update->ntimestep)
  // // 			     <<" "<<fff
  // // 			     <<", body=["<<vmax_b[0]<<","<<vmax_b[1]<<","<<vmax_b[2]<<"]"
  // // 			     <<", atom=["<<vmax_a[0]<<","<<vmax_a[1]<<","<<vmax_a[2]<<"]\n";


  // // // now body associated with each atom 

  // // // std::cout<<std::scientific;
  // // // for (int p = 0; p<comm->nprocs; p++) {
  // // //   if (p==comm->me)
  // // for (int i = 0; i < n_body(); i++) {
  // //   get_vcm(tmpv, i);
  // //   if (abs(tmpv[0])>vmax_b[0]) vmax_b[0] = abs(tmpv[0]);
  // //   if (abs(tmpv[1])>vmax_b[1]) vmax_b[1] = abs(tmpv[1]);
  // //   if (abs(tmpv[2])>vmax_b[2]) vmax_b[2] = abs(tmpv[2]);

  // //   std::vector<int> iatoms = atoms_of_body(tag(i), atom->nlocal);
  // //   for (int iatom : iatoms)  {
  // //     if (abs(v[iatom][0])>vmax_a[0]) vmax_a[0] = abs(v[iatom][0]);
  // //     if (abs(v[iatom][1])>vmax_a[1]) vmax_a[1] = abs(v[iatom][1]);
  // //     if (abs(v[iatom][2])>vmax_a[2]) vmax_a[2] = abs(v[iatom][2]);
  // //   }
  // //   // std::cout<< "\n"<<std::flush;
  // // }

  // // int check = 0;
  // // if (abs(vmax_b[0])>1e-16) check = std::max(int(vmax_a[0]/vmax_b[0]>2.0  || vmax_a[0]/vmax_b[0]<0.5),check);
  // // if (abs(vmax_b[1])>1e-16) check = std::max(int(vmax_a[1]/vmax_b[1]>2.0  || vmax_a[1]/vmax_b[1]<0.5),check);
  // // if (abs(vmax_b[2])>1e-16) check = std::max(int(vmax_a[2]/vmax_b[2]>2.0  || vmax_a[2]/vmax_b[2]<0.5),check);

  // // int c2 = 0;
  // // MPI_Reduce(&check, &c2, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

  // // if (comm->me==0 && c2) std::cout <<std::to_string(update->ntimestep)
  // // 			     <<" "<<fff
  // // 			     <<", body=["<<vmax_b[0]<<","<<vmax_b[1]<<","<<vmax_b[2]<<"]"
  // // 			     <<", atom=["<<vmax_a[0]<<","<<vmax_a[1]<<","<<vmax_a[2]<<"]\n";


  // test 1
  // // return ;
  // // if (update->ntimestep<932000) return; 
  // // if (update->ntimestep<75000) return; 
  // // if (update->ntimestep<493650) return; 
  // // if (update->ntimestep<529165) return; 
  // // if (update->ntimestep<250000) return; 
  // int dbg_tag = 19;//79;//11;//23;//315;//355;//281;//72;
  // int max_comm = 0;
  // double tmpx[3],tmpv[3];
  // tagint *image = atom->image;
  // std::cout<<std::scientific;
  // for (int p = 0; p<comm->nprocs; p++) {
  //   if (p==comm->me)
  //     for (int i = 0; i < n_body(); i++) {
  // 	// get_xcm(tmpx, i);
  //       multisphere_.x_bound(tmpx, i); // has to use this, otherwise there is a sync error from the shifts performed here and the exchange of bodies?      
  // 	get_vcm(tmpv, i);
  // 	// if ((tag(i)==dbg_tag || tag(i)==21)){
  // 	if (tag(i)==dbg_tag) {
  // 	  // std::cout <<std::to_string(update->ntimestep)<<" "<<fff<<", cpu="<< p << ", tag("<<i<<") = "<<tag(i)
  // 	  // 	    <<", x=["<<tmpx[0]<<","<<tmpx[1]<<","<<tmpx[2]<<"]";
  // 	  std::vector<int> iatoms = atoms_of_body(tag(i), atom->nlocal // + atom->nghost
  // 						  );
  // 	  double rmin = 1e13,rmax=-1e13;
  // 	  int nna = 0;
  // 	  int rma = 0;
  // 	  for (int iatom : iatoms)  {
  // 	    double xp = atom->x[iatom][0]-tmpx[0];
  // 	    double yp = atom->x[iatom][1]-tmpx[1];
  // 	    double zp = atom->x[iatom][2]-tmpx[2];
  // 	    double rr = std::sqrt(xp*xp+yp*yp+zp*zp);
  // 	    if (rr<1e-7) rma++;
  // 	    rmin = std::min(rmin,rr);
  // 	    rmax = std::max(rmax,rr);
  // 	    nna++;
  // 	    // if (xbox || ybox || zbox)imgcheck=1;
  // 	    // std::cout <<", [yat= "<<atom->x[iatom][1]<<", im="<<xbox<<"|"<<ybox<<"|"<<zbox<<"]";
  // 	  }
	  
  // 	  // if (rmax<1e-7 && nna)
  // 	  if (rma>=3) {
  // 	    print = 1;
  // 	    std::cout <<std::to_string(update->ntimestep)<<" "<<fff<<", cpu="<< p << ", tag("<<i<<") = "<<tag(i)
  // 		      <<", x=["<<tmpx[0]<<","<<tmpx[1]<<","<<tmpx[2]<<"]";	  
  // 	    std::cout <<", r=[ "<<rmin<<", "<<rmax<<" ], "<<nna;
  // 	    std::cout<< "\n"<<std::flush;
  // 	  }
  // 	  max_comm = std::max(max_comm,p);
  // 	}
  //     }
  //   MPI_Barrier(MPI_COMM_WORLD);
  // }

  // // int nlocal = atom->nlocal + atom->nghost;
  // // double **x = atom->x;
  // // double **v = atom->v;
  // // // now atoms, check that they belong to a body

  // // for (int p = 0; p<comm->nprocs; p++) {
  // //   if (p==comm->me)
  // //     for (int i = 0; i < nlocal; i++) {
  // // 	if (this->belongs_to(i)==dbg_tag // || this->belongs_to(i)>20
  // // 	    )  {
  // // 	    int xbox = (image[i] & IMGMASK) - IMGMAX;
  // // 	    int ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
  // // 	    int zbox = (image[i] >> IMG2BITS) - IMGMAX;	  

  // // 	  std::cout<<std::to_string(update->ntimestep)<<" " <<fff<<", cpu="<< p << ", belongs_to("<<i<<") = "<<  this->belongs_to(i)
  // // 		   <<", x=["<<x[i][0]<<","<<x[i][1]<<","<<x[i][2]<<"], " << "ghost="<<((i>=atom->nlocal)?1:0)
  // // 		   <<", im="<<xbox<<"|"<<ybox<<"|"<<zbox<<"\n";
  // // 	    // <<", x=["<<x[i][0]<<","<<x[i][1]<<","<<x[i][2]<<"]"
  // // 	    // <<", v=["<<v[i][0]<<","<<v[i][1]<<","<<v[i][2]<<"]\n";
  // // 	  max_comm = std::max(max_comm,p);
  // // 	}
	
  // //     }
    
  // //   MPI_Barrier(MPI_COMM_WORLD);
  // // }


  // // // method 2
  // //   ///////////////////////////////////////////////////////  
  // // return;
  // // MPI_Barrier(MPI_COMM_WORLD);
  // // if (comm->me==0 && domain->switchLEbc) std::cout <<std::to_string(update->ntimestep)<<" "<<fff<<"\n"<<std::flush;
  // // MPI_Barrier(MPI_COMM_WORLD);
  // return;
  // // static int imcheck = 0;
  // MPI_Barrier(MPI_COMM_WORLD);
  // if (update->ntimestep<1440025) return; //1149690
  // // if (update->ntimestep<2379991 || update->ntimestep>2379993) return; 
  // // if (update->ntimestep<493650) return; 
  // // if (update->ntimestep<529165) return; 
  // // if (update->ntimestep<250000) return; 
  // int dbg_tag =758;//79;//11;//23;//315;//355;//281;//72;
//   int dbg_tag =ex;
//   double tmpx[3],tmpv[3];
//   tagint *image = atom->image;
//   std::cout<<std::scientific;
//   for (int p = 0; p<comm->nprocs; p++) {
//     if (p==comm->me)
//       for (int i = 0; i < n_body(); i++) {
//   	// get_xcm(tmpx, i);
//         multisphere_.x_bound(tmpx, i); // has to use this, otherwise there is a sync error from the shifts performed here and the exchange of bodies?      
//   	get_vcm(tmpv, i);
//   	// if ((tag(i)==dbg_tag || tag(i)==21)){
//   	if (tag(i)==dbg_tag // || tag(i)>20
//   	    ) {
//   	  std::cout <<std::to_string(update->ntimestep)<<" "<<fff<<", cpu="<< p << ", tag("<<i<<") = "<<tag(i)
//   		    <<", x=["<<tmpx[0]<<","<<tmpx[1]<<","<<tmpx[2]<<"]";
//   	  // <<", x=["<<tmpx[0]<<","<<tmpx[1]<<","<<tmpx[2]<<"]"
//   	  // <<", v=["<<tmpv[0]<<","<<tmpv[1]<<","<<tmpv[2]<<"]";
//   	  std::vector<int> iatoms = atoms_of_body(tag(i), atom->nlocal + atom->nghost);
//   	  for (int iatom : iatoms)  {
//   	    // int xbox = (image[iatom] & IMGMASK) - IMGMAX;
//   	    // int ybox = (image[iatom] >> IMGBITS & IMGMASK) - IMGMAX;
//   	    // int zbox = (image[iatom] >> IMG2BITS) - IMGMAX;
//   	    // if (xbox || ybox || zbox)imgcheck=1;
//   	    std::cout <<", [yat= "<<atom->x[iatom][1];//<<", im="<<xbox<<"|"<<ybox<<"|"<<zbox<<"]";
//   	  }
//   	  std::cout<< "\n"<<std::flush;
//   	}
//       }
//     MPI_Barrier(MPI_COMM_WORLD);
//   }
//   int nlocal = atom->nlocal;// + atom->nghost;
//   double **x = atom->x;
//   double **v = atom->v;
//   // now atoms, check that they belong to a body
//   int max_comm = 0;
//   for (int p = 0; p<comm->nprocs; p++) {
//     if (p==comm->me)
//       for (int i = 0; i < nlocal; i++) {
//   	if (this->belongs_to(i)==dbg_tag // || this->belongs_to(i)>20
//   	    )  {
//   	  // int xbox = (image[i] & IMGMASK) - IMGMAX;
//   	  // int ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
//   	  // int zbox = (image[i] >> IMG2BITS) - IMGMAX;	  

//   	  std::cout<<std::to_string(update->ntimestep)<<" " <<fff<<", cpu="<< p << ", belongs_to("<<i<<") = "<<  this->belongs_to(i)
//   		   <<", x=["<<x[i][0]<<","<<x[i][1]<<","<<x[i][2]<<"], \n" ;
// 	    //<< "ghost="<<((i>=atom->nlocal)?1:0)<<", im="<<xbox<<"|"<<ybox<<"|"<<zbox<<"\n";
//   	  // <<", x=["<<x[i][0]<<","<<x[i][1]<<","<<x[i][2]<<"]"
//   	  // <<", v=["<<v[i][0]<<","<<v[i][1]<<","<<v[i][2]<<"]\n";
//   	  max_comm = std::max(max_comm,p);
//   	}
	
//       }
    
//     MPI_Barrier(MPI_COMM_WORLD);
//   }
  
//   MPI_Barrier(MPI_COMM_WORLD);
//   int mcc = 0;
//   MPI_Allreduce(&max_comm,&mcc,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
//   if (comm->me==mcc && print)std::cout<< "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n"<<std::flush;
//   MPI_Barrier(MPI_COMM_WORLD);

}

