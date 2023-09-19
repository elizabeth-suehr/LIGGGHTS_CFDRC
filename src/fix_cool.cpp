/*
// LICENSE Date: 08-26-2019
// This file is part of the CFDRC Gas-Granular Flow Solver Loci/GGFS, Copyright 2019.
// The software tool Loci/GGFS module consisting
// is being furnished to NASA Personnel under SBIR Data Rights.
// The SBIR DATA rights are asserted by CFD Research Corporation.
// Distribution C: Limited to Government Employees only.
// A release under Distribution B and A is being considered and may be done for future releases of the code.
//
// For more information, contact:
//
// Manuel P. Gale (manuel.gale@cfdrc.com), 256-726-4860.
//
// These SBIR data are furnished with SBIR rights under
// Contract No. 80NSSC18P2154. For a period of 4 years, unless
// extended in accordance with FAR 27.409(h), after acceptance
// of all items to be delivered under this contract, the
// Government will use these data for Government purposes
// only, and they shall not be disclosed outside the
// Government (including disclosure for procurement purposes)
// during such period without permission of the Contractor,
// except that, subject to the foregoing use and disclosure
// prohibitions, these data may be disclosed for use by
// support Contractors. After the protection period, the
// Government has a paid-up license to use, and to authorize
// others to use on its behalf, these data for Government
// purposes, but is relieved of all disclosure prohibitions
// and assumes no liability for unauthorized use of these data
// by third parties. This notice shall be affixed to any
// reproductions of these data, in whole or in part.
*/

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

    Andreas Aigner (JKU Linz)
    Andreas Eitzlmayr (TU Graz)

    Copyright 2009-2012 JKU Linz
    Copyright 2013-     TU Graz
------------------------------------------------------------------------- */

#include <cmath>
#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include "fix_cool.h"
#include "update.h"
#include "respa.h"
#include "atom_vec.h"
#include "atom.h"
#include "force.h"
#include "modify.h"
#include "pair.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "timer.h"
#include <vector>
#include <iostream>
#include <exception>
#include <csignal>
#include <stdexcept>
#include <limits>
#include <fenv.h>
#include "assert.h"
#include <string>
#include <sstream>
#include "fix_multisphere.h"
#include "particleToInsert_multisphere.h"
#include "fix_template_multisphere.h"
#include <random>
#include <unistd.h>
#include "fix_relax_contacts.h"  
#include <iostream>
#include <fstream>
#include <iomanip>
#include "signal_handling.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace std;

/* ---------------------------------------------------------------------- */
// #define VERBOUT
FixCOOL::FixCOOL(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg) {
  global_freq = 1;
  // Check args
  int iarg = 3;
  if (narg < iarg + 1)
    error->fix_error(FLERR, this, "Not enough arguments \n");

  if (strcmp(arg[iarg], "gtemp") == 0) {
    if (narg < iarg + 1)
      error->fix_error(FLERR, this, "Not enough arguments for 'gtemp' cooling specification \n");
    gtemp = force->numeric(FLERR, arg[iarg + 1]);
    iarg += 2;
  }

  fix_ms = static_cast<FixMultisphere *>(modify->find_fix_style("multisphere", 0));

  int multisphere_types = 0;
  if (strcmp(arg[iarg], "multisphere_types") == 0)
  {
    if (narg < iarg + 1)
      error->fix_error(FLERR, this, "Not enough arguments for 'multisphere_types' cooling specification \n");
    multisphere_types = force->inumeric(FLERR, arg[iarg + 1]);
    iarg += 2;
  }
  
//   for (int i = 0; i < multisphere_types; i++){
//     fix_ms_templates.push_back(static_cast<FixTemplateMultisphere *>(modify->find_fix_style("particletemplate/multisphere", i)));
//     if (!fix_ms_templates[i]) {
//       std::cout << "incorrect particletemplates/multisphere" << std::endl;
//       exit(0);
//     }
//   }

  lattice = 0;
  if (strcmp(arg[iarg], "lattice") == 0)
  {
    if (narg < iarg + 1)
      error->fix_error(FLERR, this, "Not enough arguments for 'lattice' cooling specification \n");
    lattice = force->numeric(FLERR, arg[iarg + 1]);
    iarg += 2;
  }


  // EMS COMMENT: Might be able to remove?
  atom->add_callback(0);

  if (comm->me == 0)
    cout << "FixCOOL::init ... granular_temperature=" << gtemp
	 << ", body_types=" << multisphere_types
	 << ", lattice=" << lattice << endl;
}

/* ---------------------------------------------------------------------- */

FixCOOL::~FixCOOL()
{
  atom->delete_callback(id, 0);
}

/* ---------------------------------------------------------------------- */

int FixCOOL::setmask()
{
  int mask = 0;
  // mask |= END_OF_STEP;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

class loc {
public:
  double x,y,z;
};

void FixCOOL::init() {

  SignalHandler::request_exit=false;
  if (comm->me==0) std::cout << "Fix::init()" << std::endl;
  Fix::init();
  if (comm->me==0) std::cout << "FixCOOL::init()" << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);
  atom->lebc_highest_atom_id = atom->natoms;

  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(-1,1);
  
  
 {
    double v = sqrt(3.*gtemp);
    if (fix_ms) {
      
      // EMS COMMENT: sets up access to fix_ms in domain and comm
      fix_ms->lebc_highest_body_id = fix_ms->tag_max_body();
      for (int i = 0; i < fix_ms->n_body(); i++)
	{
	  // fix_ms->data().x_bound(x_cm, i);

	  // gtemp = 1/3*u.u, umag =sqrt(gtemp)
      	  double randN = distribution(generator);
      	  double randN_x = ((distribution(generator)>0.0)?1.0:-1.0)*randN;
      	  double randN_y = ((distribution(generator)>0.0)?1.0:-1.0)*randN;
      	  double randN_z = ((distribution(generator)>0.0)?1.0:-1.0)*randN;
      	  double v_cm[3] = {randN*v, randN*v, randN*v};
      	  // double v_cm[3] = {0, 0.0, 0};

	  // double v_x = domain->ssr * (x_cm[1] - ylength * 0.5);

	  // double randN = distribution(generator);

	  // double v_cm[3] = {v_x, 0.3*randN*v_x, 0.3*randN*v_x};
	  // // double v_cm[3] = {0, -1.0, 0};

	  fix_ms->data().set_v_body(i, v_cm);
	  double cm[3] = {0};
	  fix_ms->data().reset_forces(0);
	  fix_ms->data().set_angmom_via_omega_body(i, cm);
	}
      for (int i = 0; i < atom->nlocal; i++)
	{
	  if (fix_ms->belongs_to(i) == -1)
	    {
	      double randN = distribution(generator);
	      double randN_x = ((distribution(generator)>0.0)?1.0:-1.0)*randN;
	      double randN_y = ((distribution(generator)>0.0)?1.0:-1.0)*randN;
	      double randN_z = ((distribution(generator)>0.0)?1.0:-1.0)*randN;
	      double v_cm[3] = {randN*v, randN*v, randN*v};

	      atom->v[i][0] =  v_cm[0];
	      atom->v[i][1] =  v_cm[1];
	      atom->v[i][2] =  v_cm[2];
	    }
	}
    }
    else
      {
	for (int i = 0; i < atom->nlocal; i++) {
	  // gtemp = 1/3*u.u, umag =sqrt(gtemp)
	  double randN = distribution(generator);
	  double randN_x = ((distribution(generator)>0.0)?1.0:-1.0)*randN;
	  double randN_y = ((distribution(generator)>0.0)?1.0:-1.0)*randN;
	  double randN_z = ((distribution(generator)>0.0)?1.0:-1.0)*randN;
	  double v_cm[3] = {randN*v, randN*v, randN*v};

	  atom->v[i][0] =  v_cm[0];
	  atom->v[i][1] =  v_cm[1];
	  atom->v[i][2] =  v_cm[2];


	}
      }
    
  }

  if (fix_ms) { // needed
    fix_ms->set_xv();
    fix_ms->out_rev_comm_x_v_omega();
    
    fix_ms->add_body_finalize();
    // // clear map and since we are deleting all of the ghost particles
    // atom->nghost = 0;
    // atom->map_init();
    // atom->map_set();  
    // // fix_ms->add_body_finalize();
    
    // comm->borders();      
    // // // neighbor->build();
    
    // fix_ms->comm_after_mod();
  }

}




void FixCOOL::end_of_step() {

  // if (update->ntimestep % 1000 == 0) {
  //   double cke = check_ke(); 
  //   if (comm->me == 0) std::cout << "step: "<< update->ntimestep
  // 				 <<", natom="<<atom->natoms
  // 				 <<", nbody="<<(fix_ms?fix_ms->n_body_all():0)
  // 				 <<", ke="<<cke
  // 				 << std::endl;
  // }

}

void FixCOOL::post_force(int dummy) {

  int *mask = atom->mask;
  int *itype = atom->type;
  double *density = atom->density;
  double *radius = atom->radius;
  // double *p = atom->p;
  int nlocal = atom->nlocal;
  int nglobal = atom->natoms;

  double **v = atom->v;

  double local_sum_v[3] = {0, 0, 0};
  double domain_volume = domain->xprd * domain->yprd * domain->zprd;

  double mass_sum = 0.0;
  if (atom->rmass)
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit && (atom->tag[i] <= atom->lebc_highest_atom_id && ((fix_ms)?(fix_ms->belongs_to(i) == -1):1))) {
	double massp = atom->rmass[i];
	local_sum_v[0] += v[i][0] * massp;
	local_sum_v[1] += v[i][1] * massp;
	local_sum_v[2] += v[i][2] * massp;
	mass_sum += massp;
      }
    }
  else 
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit ) {
	double massp = atom->mass[itype[i]];
	local_sum_v[0] += v[i][0] * massp;
	local_sum_v[1] += v[i][1] * massp;
	local_sum_v[2] += v[i][2] * massp;
	mass_sum += massp;
      }
    }
      
  
  // double cm[3] = {0};
  if (fix_ms)
    for (int i = 0; i < fix_ms->n_body(); i++){
      if (fix_ms->data().tag(i) != -1){
	double velocity[3] = {0.0};
	// fix_ms->get_xcm(cm,i);
	// if (cm[1] < lo[1] || cm[1] > hi[1]) continue;
        fix_ms->data().vcm(velocity,i);
	double massB = fix_ms->data().mass(i);
        local_sum_v[0] += velocity[0] * massB;
        local_sum_v[1] += velocity[1] * massB;
        local_sum_v[2] += velocity[2] * massB;
        mass_sum += massB;
      }
    }
  
  double global_sum_v[3];
  MPI_Allreduce(&local_sum_v[0], &global_sum_v[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  double global_sum_mass;
  MPI_Allreduce(&mass_sum, &global_sum_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  double global_avg_v[3];
  global_avg_v[0] = global_sum_v[0] / global_sum_mass;
  global_avg_v[1] = global_sum_v[1] / global_sum_mass;
  global_avg_v[2] = global_sum_v[2] / global_sum_mass;
  
  double k_stress_tensor_sum[6] = {0, 0, 0, 0, 0, 0};
  
  if (atom->rmass) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit && (atom->tag[i] <= atom->lebc_highest_atom_id && ((fix_ms)?(fix_ms->belongs_to(i) == -1):1))) {
	double mss = atom->rmass[i];
	k_stress_tensor_sum[0] += mss * (v[i][0] - global_avg_v[0]) * (v[i][0] - global_avg_v[0]); // xx
	k_stress_tensor_sum[1] += mss * (v[i][1] - global_avg_v[1]) * (v[i][1] - global_avg_v[1]); // yy
	k_stress_tensor_sum[2] += mss * (v[i][2] - global_avg_v[2]) * (v[i][2] - global_avg_v[2]); // zz
	k_stress_tensor_sum[3] += mss * (v[i][0] - global_avg_v[0]) * (v[i][1] - global_avg_v[1]); // xy
	k_stress_tensor_sum[4] += mss * (v[i][0] - global_avg_v[0]) * (v[i][2] - global_avg_v[2]); // xz
	k_stress_tensor_sum[5] += mss * (v[i][1] - global_avg_v[1]) * (v[i][2] - global_avg_v[2]); // yz
      }
    }
  } else {      
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	double mss = atom->mass[itype[i]];
	k_stress_tensor_sum[0] += mss * (v[i][0] - global_avg_v[0]) * (v[i][0] - global_avg_v[0]); // xx
	k_stress_tensor_sum[1] += mss * (v[i][1] - global_avg_v[1]) * (v[i][1] - global_avg_v[1]); // yy
	k_stress_tensor_sum[2] += mss * (v[i][2] - global_avg_v[2]) * (v[i][2] - global_avg_v[2]); // zz
	k_stress_tensor_sum[3] += mss * (v[i][0] - global_avg_v[0]) * (v[i][1] - global_avg_v[1]); // xy
	k_stress_tensor_sum[4] += mss * (v[i][0] - global_avg_v[0]) * (v[i][2] - global_avg_v[2]); // xz
	k_stress_tensor_sum[5] += mss * (v[i][1] - global_avg_v[1]) * (v[i][2] - global_avg_v[2]); // yz	  
      }      
  }
  
  if (fix_ms)
    for (int i = 0; i < fix_ms->n_body(); i++){
      if (fix_ms->data().tag(i) != -1) {
	double velocity[3] = {0.0};
        fix_ms->data().vcm(velocity,i);
    	double masss = fix_ms->data().mass(i);
        k_stress_tensor_sum[0] += masss * (velocity[0] - global_avg_v[0]) * (velocity[0] - global_avg_v[0]); // xx
        k_stress_tensor_sum[1] += masss * (velocity[1] - global_avg_v[1]) * (velocity[1] - global_avg_v[1]); // yy
        k_stress_tensor_sum[2] += masss * (velocity[2] - global_avg_v[2]) * (velocity[2] - global_avg_v[2]); // zz
        k_stress_tensor_sum[3] += masss * (velocity[0] - global_avg_v[0]) * (velocity[1] - global_avg_v[1]); // xy
        k_stress_tensor_sum[4] += masss * (velocity[0] - global_avg_v[0]) * (velocity[2] - global_avg_v[2]); // xz
        k_stress_tensor_sum[5] += masss * (velocity[1] - global_avg_v[1]) * (velocity[2] - global_avg_v[2]); // yz
      }
    }

    double global_k_stress_tensor_sum[6];
    MPI_Allreduce(&k_stress_tensor_sum[0], &global_k_stress_tensor_sum[0], 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double global_avg_k_stress_tensor[6];

    for (int i = 0; i < 6; i++)
    {
      global_avg_k_stress_tensor[i] = global_k_stress_tensor_sum[i] / domain_volume;
    }

    for (int i = 0; i < 6; i++)
    {
      global_time_avg_k_stress_tensor[i] = (global_time_avg_k_stress_tensor[i] * averaging_count + global_avg_k_stress_tensor[i]) / (averaging_count + 1);
    }
    averaging_count++;

    int NN = 10000;
    if (averaging_count == NN || update->ntimestep==1)
      {
	averaging_count = 0;

      std::ostringstream oss;
      double Pk = global_time_avg_k_stress_tensor[1];
      double Pc = 0;
      double sigk = global_time_avg_k_stress_tensor[3];
      double sigc = 0;
      double time_out = update->dt*update->ntimestep;
      double P_k = (1./3.*(global_time_avg_k_stress_tensor[0] + 
			   global_time_avg_k_stress_tensor[1] + 
			   global_time_avg_k_stress_tensor[2]) ); // in this case = sum(rho_i*alpha_i) * theta
      // sum(rho_i*alpha_i) = sum(rho_i*volp_i)/V_domain = sum(mass_i)/V_domain
      double theta_sum = P_k/global_sum_mass * domain_volume;

      if (first_step) {
	theta0 = theta_sum;
	time0 = time_out;
	first_step=0;
      }

      ///////////////////////// test to see if we can predict when to stop the simulation
      double tmp_dTdt[10] = {0};
      double tmp_dTdtl[10] = {0};
      for (int i=0; i<10; i++) {
	tmp_dTdt[i] = dTdt[i];
	tmp_dTdtl[i] = dTdtl[i];
      }

      for (int i=1; i<10; i++) {
        dTdt[i] = tmp_dTdt[i-1];
	dTdtl[i] = tmp_dTdtl[i-1];
      }
#define fun(x) log10(x)
      // derivatives
      double delt = NN*update->dt;
      dTdtl[0]   = (1.5*fun(theta_sum/theta0) - 2.0*thefl[0] + 0.5*thefl[1])/(fun(1.0+delt/(update->dt*(update->ntimestep-1))));
      dTdt[0]    = (1.5*(theta_sum) - 2.0*thef[0] + 0.5*thef[1])/delt;
      // d2Tdt2[0] = (2.0*fun(theta_sum/theta0) - 5.0*thef[0] + 4.0*thef[1] - thef[2])/delt/delt;

      double tmp_thef[10] = {0};
      double tmp_thefl[10] = {0};
      for (int i=0; i<3; i++) {
	tmp_thefl[i] = thefl[i];
	tmp_thef[i] = thef[i];
      }

      for (int i=1; i<3; i++) {
        thefl[i] = tmp_thefl[i-1];
        thef[i] = tmp_thef[i-1];
      }
      thefl[0] = fun(theta_sum/theta0);
      thef[0] = theta_sum;

      double dtdt_mean = 0;
      for (int i=0; i<10; i++)
	dtdt_mean += dTdtl[i];
      dtdt_mean /= 10.0;

      if (comm->me==0)oss <<std::scientific<<std::setprecision(3)<< "theta: "<< time_out
	  <<" " << theta_sum
	  <<" " << P_k 
	  <<" " << dTdt[0]
	  <<" " << dTdtl[0]
	  <<" " << dtdt_mean << "\n";

      std::string text = oss.str();

      char *char_arr;
      char_arr = &text[0];
      // std::cout << char_arr;
      if (comm->me==0) {
        cout << char_arr;
        fprintf(logfile, "%.6e %.6e %.6e %.6e %.6e %.6e\n", time_out,theta_sum,P_k,dTdt[0],dTdtl[0],dtdt_mean);
        fflush(logfile);
      }

      for (int i = 0; i < 6; i++)
      {
        global_time_avg_k_stress_tensor[i] = 0;
      }

    }

}

