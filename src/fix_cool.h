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

#ifdef FIX_CLASS

FixStyle(cool, FixCOOL)

#else

#ifndef LMP_FIX_COOL_H
#define LMP_FIX_COOL_H

#include "fix.h"
#include "domain.h"
#include "particleToInsert.h"

namespace LAMMPS_NS
{

  class FixCOOL : public Fix
  {
  public:
    FixCOOL(class LAMMPS *, int, char **);
    ~FixCOOL();
    int setmask();
    void init();
    void post_force(int dummy);
    void end_of_step();
    /* double check_ke(); */

  private:
    double gtemp=0.0;
    int lattice=0;
    double global_time_avg_k_stress_tensor[6] = {0};
    double thef[3] = {0};
    double thefl[3] = {0};
    double dTdt[10] = {0};
    double dTdtl[10] = {0};
    /* double d2Tdt2[10] = {0}; */
    int first_step=1;
    double theta0=1.0;
    double time0=1.0;
    int averaging_count = 0;
    /* ----------------------------------------------------------------------*/
    class FixMultisphere *fix_ms;

  };

}

#endif
#endif
