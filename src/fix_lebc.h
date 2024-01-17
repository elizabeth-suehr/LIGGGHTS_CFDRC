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

#ifdef FIX_CLASS

FixStyle(lebc, FixLEBC)

#else

#ifndef LMP_FIX_LEBC_H
#define LMP_FIX_LEBC_H

#include "fix.h"
#include "fix_multisphere.h"

namespace LAMMPS_NS
{

    class FixLEBC : public Fix
    {
    public:
        FixLEBC(class LAMMPS *, int, char **);
        ~FixLEBC();
        int setmask();
        void init();
        void setup(int);
        void initial_integrate(int vflag);
        void post_integrate();
        void post_force(int);
        void print_body_data();

        int pack_border(int n, int *list, double *buf, int *pbc);
        int unpack_border(int n, int first, double *buf);
        void pre_neighbor();
        void forward_comm();
        void reverse_comm();

        int pack_comm(int, int *, double *, int, int *);
        int pack_reverse_comm(int, int, double *);
        void unpack_comm(int n, int first, double *buf);
        void unpack_comm_f_test(int n, int first, double *buf);

        int pack_comm_f_test(int, int *, double *, int, int *);
        int pack_reverse_comm_f_test(int, int, double *);

    private:
        double ssr = 0.0;
        double kinetic_stress_tensor[6] = {0.0};
        double collision_stress_tensor[6] = {0.0};
        int ave_count = 0;
        int ave_count_reset = 50000;
        bool isDimensional = true;

        FixMultisphere *fix_ms;
        int init_total_bodies = 0;

        double gtemp_distribution = 0.000001;

        string body_data_name;
        int save_count = 0;
        int save_count_reset = 5000;
    };

}

#endif
#endif
