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

#ifndef LMP_MULTISPHERE_PARALLEL_H
#define LMP_MULTISPHERE_PARALLEL_H

#define LMP_MULTISPHERE_PARALLEL_FLAG_

#include "multisphere.h"
#include "comm.h"

namespace LAMMPS_NS
{

  class MultisphereParallel : public Multisphere
  {

  public:
    MultisphereParallel(LAMMPS *lmp);
    ~MultisphereParallel();

    void exchange();

    void writeRestart(FILE *fp);
    void restart(double *list);

  private:
    int pack_exchange_rigid(int i, double *buf);
    int pack_exchange_rigid_lebc(double *buf, double *buf_from);
    int unpack_exchange_rigid(double *buf);

    void grow_send(int, int);
    void grow_recv(int);
    void migrate_bodies();
    int coord2proc(double *x, int &igx, int &igy, int &igz);
    int create_body(int n, int *sizes, int *proclist);

     // plan for irregular communication of atoms
    // no params refer to atoms copied to self

    

    

    // current size of send/recv buffer
    // send buffer and recv buffer for all comm
    int maxsend_, maxrecv_;
    double *buf_send_;
    double *buf_recv_;



    struct PlanBody {
      int nsend;                 // # of messages to send
      int nrecv;                 // # of messages to recv
      int sendmax;               // # of doubles in largest send message
      int *proc_send;            // procs to send to
      int *length_send;          // # of doubles to send to each proc
      int *num_send;             // # of atoms to send to each proc
      int *index_send;           // list of which atoms to send to each proc
      int *offset_send;          // where each atom starts in send buffer
      int *proc_recv;            // procs to recv from
      int *length_recv;          // # of doubles to recv from each proc
      MPI_Request *request;      // MPI requests for posted recvs
      MPI_Status *status;        // MPI statuses for WaitAll
    };

    PlanBody *aplan;

    void exchange_body(double *sendbuf, int *sizes, double *recvbuf);

    void destroy_body();
  };

// *************************************
#include "multisphere_parallel_I.h"
  // *************************************
}

#endif
