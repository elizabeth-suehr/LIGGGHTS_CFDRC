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

#define DELTA 10000

#include "multisphere_parallel.h"
#include "atom.h"
#include "atom_vec.h"
#include "vector_liggghts.h"
#include "domain.h"
#include "memory.h"

#include <iostream>

#define BUFFACTOR 1.5
#define BUFMIN 1000
#define BUFEXTRA 1000

/* ----------------------------------------------------------------------
   constructor / destructor
------------------------------------------------------------------------- */

MultisphereParallel::MultisphereParallel(LAMMPS *lmp) : Multisphere(lmp),

                                                        // initialize comm buffers & exchange memory
                                                        maxsend_(BUFMIN),
                                                        maxrecv_(BUFMIN),
                                                        buf_send_((double *)memory->smalloc((maxsend_ + BUFEXTRA) * sizeof(double), "frm:buf_send_")),
                                                        buf_recv_((double *)memory->smalloc((maxsend_ + BUFEXTRA) * sizeof(double), "frm:buf_send_"))
{

  aplan = NULL;
}

MultisphereParallel::~MultisphereParallel()
{
  memory->sfree(buf_send_);
  memory->sfree(buf_recv_);
}

/* ----------------------------------------------------------------------
   realloc the size of the send buffer as needed with BUFFACTOR & BUFEXTRA
   if flag = 1, realloc
   if flag = 0, don't need to realloc with copy, just free/malloc
------------------------------------------------------------------------- */

void MultisphereParallel::grow_send(int n, int flag)
{
  maxsend_ = static_cast<int>(BUFFACTOR * n);
  if (flag)
    buf_send_ = (double *)memory->srealloc(buf_send_, (maxsend_ + BUFEXTRA) * sizeof(double), "comm:buf_send_");
  else
  {
    memory->sfree(buf_send_);
    buf_send_ = (double *)memory->smalloc((maxsend_ + BUFEXTRA) * sizeof(double), "comm:buf_send_");
  }
}

/* ----------------------------------------------------------------------
   free/malloc the size of the recv buffer as needed with BUFFACTOR
------------------------------------------------------------------------- */

void MultisphereParallel::grow_recv(int n)
{
  maxrecv_ = static_cast<int>(BUFFACTOR * n);
  memory->sfree(buf_recv_);
  buf_recv_ = (double *)memory->smalloc(maxrecv_ * sizeof(double), "comm:buf_recv_");
}

/* ----------------------------------------------------------------------
   exchange bodies with neighbor procs
------------------------------------------------------------------------- */

void MultisphereParallel::exchange()
{

  if(domain->lebc)
  {
    migrate_bodies();
    return;
  }
 
  int i, m, nsend, nrecv, nrecv1, nrecv2;
  double lo, hi, value;
  double x[3];
  double *sublo, *subhi, *buf;
  MPI_Request request;
  MPI_Status status;

  // subbox bounds for orthogonal
  // triclinic not implemented

  sublo = domain->sublo;
  subhi = domain->subhi;

  std::vector<int> extra;

  bool is_top_processor = false;
  bool is_bottom_processor = false;

  if (domain->lebc)
  {
    if (domain->sublo[1] == domain->boundary[1][0])
    {
      is_bottom_processor = true;
    }
    if (domain->subhi[1] == domain->yprd)
    {
      is_top_processor = true;
    }
  }

  // loop over dimensions

  int dim_order[3] = {1, 0, 2};

  for (int dim_o = 0; dim_o < 3; dim_o++)
  {
    int dim = dim_order[dim_o];

    // fill buffer with atoms leaving my box, using < and >=
    // when atom is deleted, fill it in with last atom

    lo = sublo[dim];
    hi = subhi[dim];
    i = nsend = 0;

    while (i < nbody_)
    {

      MathExtraLiggghts::local_coosys_to_cartesian(x, xcm_to_xbound_(i), ex_space_(i), ey_space_(i), ez_space_(i));
      vectorAdd3D(xcm_(i), x, x);

      if (x[dim] < lo || x[dim] >= hi)
      {
        if (nsend > maxsend_)
          grow_send(nsend, 1);
        nsend += pack_exchange_rigid(i, &buf_send_[nsend]);

        remove_body(i);
      }
      else
        i++;
    }

    // send/recv atoms in both directions
    // if 1 proc in dimension, no send/recv, set recv buf to send buf
    // if 2 procs in dimension, single send/recv
    // if more than 2 procs in dimension, send/recv to both neighbors

    int procneigh[3][2];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 2; j++)
        procneigh[i][j] = comm->procneigh[i][j];

    int *procgrid = comm->procgrid;

    if (procgrid[dim] == 1)
    {
      nrecv = nsend;
      buf = buf_send_;
    }
    else
    {
      MPI_Sendrecv(&nsend, 1, MPI_INT, procneigh[dim][0], 0, &nrecv1, 1, MPI_INT, procneigh[dim][1], 0, world, &status);
      nrecv = nrecv1;
      if (procgrid[dim] > 2)
      {
        MPI_Sendrecv(&nsend, 1, MPI_INT, procneigh[dim][1], 0, &nrecv2, 1, MPI_INT, procneigh[dim][0], 0, world, &status);
        nrecv += nrecv2;
      }

      if (nrecv > maxrecv_)
        grow_recv(nrecv);

      MPI_Irecv(buf_recv_, nrecv1, MPI_DOUBLE, procneigh[dim][1], 0, world, &request);
      MPI_Send(buf_send_, nsend, MPI_DOUBLE, procneigh[dim][0], 0, world);
      MPI_Wait(&request, &status);

      if (procgrid[dim] > 2)
      {
        MPI_Irecv(&buf_recv_[nrecv1], nrecv2, MPI_DOUBLE, procneigh[dim][0], 0, world, &request);
        MPI_Send(buf_send_, nsend, MPI_DOUBLE, procneigh[dim][1], 0, world);
        MPI_Wait(&request, &status);
      }

      buf = buf_recv_;
    }

    // check incoming atoms to see if they are in my box
    // if so, add to my list

    m = 0;
    extra.clear();
    while (m < nrecv)
    {
      value = buf[m + dim + 1];

      if (value >= lo && value < hi)
        m += unpack_exchange_rigid(&buf[m]);
      else
      {
        if (domain->lebc && dim == 0 && (is_bottom_processor || is_top_processor))
          extra.push_back(m);
        m += static_cast<int>(buf[m]);
      }
    }

    if (domain->lebc && (is_bottom_processor || is_top_processor) && dim == 0)
    {

      int extra_exchange = procgrid[dim]/2 -1; // Thanks Ben Armstrong !

      for (int z = 0; z < extra_exchange; z++)
      {
        lo = sublo[dim];
        hi = subhi[dim];
        nsend = 0;

        for (int ex : extra)
        {
          if (nsend > maxsend_)
            grow_send(nsend, 1);
          nsend += pack_exchange_rigid_lebc(&buf_send_[nsend], &buf[ex]);
        }
        extra.clear();

        MPI_Sendrecv(&nsend, 1, MPI_INT, procneigh[dim][0], 0, &nrecv1, 1, MPI_INT, procneigh[dim][1], 0, world, &status);
        nrecv = nrecv1;
        if (z < extra_exchange - 1 || procgrid[dim] % 2 == 1)
        {
          MPI_Sendrecv(&nsend, 1, MPI_INT, procneigh[dim][1], 0, &nrecv2, 1, MPI_INT, procneigh[dim][0], 0, world, &status);
          nrecv += nrecv2;
        }

        if (nrecv > maxrecv_)
          grow_recv(nrecv);

        MPI_Irecv(buf_recv_, nrecv1, MPI_DOUBLE, procneigh[dim][1], 0, world, &request);
        MPI_Send(buf_send_, nsend, MPI_DOUBLE, procneigh[dim][0], 0, world);
        MPI_Wait(&request, &status);

        if (z < extra_exchange - 1 || procgrid[dim] % 2 == 1)
        {
          MPI_Irecv(&buf_recv_[nrecv1], nrecv2, MPI_DOUBLE, procneigh[dim][0], 0, world, &request);
          MPI_Send(buf_send_, nsend, MPI_DOUBLE, procneigh[dim][1], 0, world);
          MPI_Wait(&request, &status);
        }

        buf = buf_recv_;

        m = 0;

        while (m < nrecv)
        {

          value = buf[m + dim + 1];

          if (value >= lo && value < hi)
          {
            
            m += unpack_exchange_rigid(&buf[m]);
          }

          else
          {
            
            extra.push_back(m);
            m += static_cast<int>(buf[m]);
          }
        }
      }
    }
  }

  // calc_nbody_all();
  // if (total_bodies < n_body_all())
  // {
  //   std::cout << "from multsphere_parallel" << std::fflush;
  // }
}


void MultisphereParallel::migrate_bodies()
{

  // subbox bounds for orthogonal or triclinic box
  // other comm/domain data used by coord2proc()

  // subbox bounds for orthogonal
  // triclinic not implemented

  double *sublo = domain->sublo;
  double *subhi = domain->subhi;

  
  
  double *boxlo = domain->boxlo;
  double *prd = domain->prd;

  // loop over atoms, flag any that are not in my sub-box
  // fill buffer with atoms leaving my box, using < and >=
  // assign which proc it belongs to via coord2proc()
  // if coord2proc() returns me, due to round-off
  //   in triclinic x2lamda(), then keep atom and don't send
  // when atom is deleted, fill it in with last atom

  
  double x[3];

  int nsend = 0;
  int nsendatom = 0;
  int *sizes = new int[nbody_];
  int *proclist = new int[nbody_];
  int igx,igy,igz;

  int i = 0;
  while (i < nbody_) {
    MathExtraLiggghts::local_coosys_to_cartesian(x, xcm_to_xbound_(i), ex_space_(i), ey_space_(i), ez_space_(i));
    vectorAdd3D(xcm_(i), x, x);

    if (x[0] < sublo[0] || x[0] >= subhi[0] ||
        x[1] < sublo[1] || x[1] >= subhi[1] ||
        x[2] < sublo[2] || x[2] >= subhi[2]) {
      proclist[nsendatom] = coord2proc(x,igx,igy,igz);
      if (proclist[nsendatom] != comm->me) {
        if (nsend > maxsend_) grow_send(nsend,1);
        sizes[nsendatom] = pack_exchange_rigid(i,&buf_send_[nsend]);
        nsend += sizes[nsendatom];
        nsendatom++;
        remove_body(i);
      } else i++;
    } else i++;
  }

  // create irregular communication plan, perform comm, destroy plan
  // returned nrecv = size of buffer needed for incoming atoms

  int nrecv = create_body(nsendatom,sizes,proclist);
  if (nrecv > maxrecv_) grow_recv(nrecv);
  exchange_body(buf_send_,sizes,buf_recv_);
  destroy_body();

  delete [] sizes;
  delete [] proclist;

  // add received atoms to my list

  int m = 0;
  while (m < nrecv) m += unpack_exchange_rigid(&buf_recv_[m]);

}


int MultisphereParallel::coord2proc(double *x, int &igx, int &igy, int &igz)
{
  
  igx = static_cast<int> (comm->procgrid[0] * (x[0]-domain->boxlo[0]) / domain->prd[0]);
  igy = static_cast<int> (comm->procgrid[1] * (x[1]-domain->boxlo[1]) / domain->prd[1]);
  igz = static_cast<int> (comm->procgrid[2] * (x[2]-domain->boxlo[2]) / domain->prd[2]);
   
  // if (igx < 0) igx = 0;
  // if (igx >= comm->procgrid[0]) igx = comm->procgrid[0] - 1;
  // if (igy < 0) igy = 0;
  // if (igy >= comm->procgrid[1]) igy = comm->procgrid[1] - 1;
  // if (igz < 0) igz = 0;
  // if (igz >= comm->procgrid[2]) igz = comm->procgrid[2] - 1;

  return comm->grid2proc[igx][igy][igz];
}


int MultisphereParallel::create_body(int n, int *sizes, int *proclist)
{
  int i;

  // allocate plan and work vectors

  if (aplan) destroy_body();
  aplan = (PlanBody *) memory->smalloc(sizeof(PlanBody),"multisphereparallel:aplan");
  int *list = new int[comm->nprocs];
  int *count = new int[comm->nprocs];

  // nrecv = # of messages I receive

  for (i = 0; i < comm->nprocs; i++) {
    list[i] = 0;
    count[i] = 1;
  }
  for (i = 0; i < n; i++) list[proclist[i]] = 1;

  int nrecv;
  MPI_Reduce_scatter(list,&nrecv,count,MPI_INT,MPI_SUM,world);

  // allocate receive arrays

  int *proc_recv = new int[nrecv];
  int *length_recv = new int[nrecv];
  MPI_Request *request = new MPI_Request[nrecv];
  MPI_Status *status = new MPI_Status[nrecv];

  // nsend = # of messages I send

  for (i = 0; i < comm->nprocs; i++) list[i] = 0;
  for (i = 0; i < n; i++) list[proclist[i]] += sizes[i];

  int nsend = 0;
  for (i = 0; i < comm->nprocs; i++)
    if (list[i]) nsend++;

  // allocate send arrays

  int *proc_send = new int[nsend];
  int *length_send = new int[nsend];
  int *num_send = new int[nsend];
  int *index_send = new int[n];
  int *offset_send = new int[n];

  // list still stores size of message for procs I send to
  // proc_send = procs I send to
  // length_send = total size of message I send to each proc
  // to balance pattern of send messages:
  //   each proc begins with iproc > me, continues until iproc = me
  // reset list to store which send message each proc corresponds to

  int iproc = comm->me;
  int isend = 0;
  for (i = 0; i < comm->nprocs; i++) {
    iproc++;
    if (iproc == comm->nprocs) iproc = 0;
    if (list[iproc] > 0) {
      proc_send[isend] = iproc;
      length_send[isend] = list[iproc];
      list[iproc] = isend;
      isend++;
    }
  }

  // num_send = # of atoms I send to each proc

  for (i = 0; i < nsend; i++) num_send[i] = 0;
  for (i = 0; i < n; i++) {
    isend = list[proclist[i]];
    num_send[isend]++;
  }

  // count = offsets into index_send for each proc I send to
  // index_send = list of which atoms to send to each proc
  //   1st N1 values are atom indices for 1st proc,
  //   next N2 values are atom indices for 2nd proc, etc
  // offset_send = where each atom starts in send buffer

  count[0] = 0;
  for (i = 1; i < nsend; i++) count[i] = count[i-1] + num_send[i-1];

  for (i = 0; i < n; i++) {
    isend = list[proclist[i]];
    index_send[count[isend]++] = i;
    if (i) offset_send[i] = offset_send[i-1] + sizes[i-1];
    else offset_send[i] = 0;
  }

  // tell receivers how much data I send
  // sendmax = largest # of doubles I send in a single message

  int sendmax = 0;
  for (i = 0; i < nsend; i++) {
    MPI_Send(&length_send[i],1,MPI_INT,proc_send[i],0,world);
    sendmax = MAX(sendmax,length_send[i]);
  }

  // receive incoming messages
  // proc_recv = procs I recv from
  // length_recv = total size of message each proc sends me
  // nrecvsize = total size of data I recv

  int nrecvsize = 0;
  for (i = 0; i < nrecv; i++) {
    MPI_Recv(&length_recv[i],1,MPI_INT,MPI_ANY_SOURCE,0,world,status);
    proc_recv[i] = status->MPI_SOURCE;
    nrecvsize += length_recv[i];
  }

  // barrier to insure all MPI_ANY_SOURCE messages are received
  // else another proc could proceed to exchange_body() and send to me

  MPI_Barrier(world);

  // free work vectors

  delete [] count;
  delete [] list;

  // initialize plan

  aplan->nsend = nsend;
  aplan->nrecv = nrecv;
  aplan->sendmax = sendmax;

  aplan->proc_send = proc_send;
  aplan->length_send = length_send;
  aplan->num_send = num_send;
  aplan->index_send = index_send;
  aplan->offset_send = offset_send;
  aplan->proc_recv = proc_recv;
  aplan->length_recv = length_recv;

  aplan->request = request;
  aplan->status = status;

  return nrecvsize;
}


void MultisphereParallel::exchange_body(double *sendbuf, int *sizes, double *recvbuf)
{
  int i,m,n,offset,num_send;

  // post all receives

  offset = 0;
  for (int irecv = 0; irecv < aplan->nrecv; irecv++) {
    MPI_Irecv(&recvbuf[offset],aplan->length_recv[irecv],MPI_DOUBLE,
              aplan->proc_recv[irecv],0,world,&aplan->request[irecv]);
    offset += aplan->length_recv[irecv];
  }

  // allocate buf for largest send

  double *buf;
  memory->create(buf,aplan->sendmax,"multisphereparallel:buf");

  // send each message
  // pack buf with list of atoms
  // m = index of atom in sendbuf

  int *index_send = aplan->index_send;
  int nsend = aplan->nsend;
  n = 0;

  for (int isend = 0; isend < nsend; isend++) {
    offset = 0;
    num_send = aplan->num_send[isend];
    for (i = 0; i < num_send; i++) {
      m = index_send[n++];
      memcpy(&buf[offset],&sendbuf[aplan->offset_send[m]],
             sizes[m]*sizeof(double));
      offset += sizes[m];
    }
    MPI_Send(buf,aplan->length_send[isend],MPI_DOUBLE,
             aplan->proc_send[isend],0,world);
  }

  // free temporary send buffer

  memory->destroy(buf);

  // wait on all incoming messages

  if (aplan->nrecv) MPI_Waitall(aplan->nrecv,aplan->request,aplan->status);
}


void MultisphereParallel::destroy_body()
{
  delete [] aplan->proc_send;
  delete [] aplan->length_send;
  delete [] aplan->num_send;
  delete [] aplan->index_send;
  delete [] aplan->offset_send;
  delete [] aplan->proc_recv;
  delete [] aplan->length_recv;
  delete [] aplan->request;
  delete [] aplan->status;
  memory->sfree(aplan);
  aplan = NULL;
}




/* ----------------------------------------------------------------------
   restart functionality - write all required data into restart buffer
   executed on all processes, but only proc 0 writes into writebuf
------------------------------------------------------------------------- */

void MultisphereParallel::writeRestart(FILE *fp)
{
  double *sendbuf = 0, *recvbuf = 0;
  double xbnd[3];
  bool dummy = false;
  double nba = static_cast<double>(n_body_all());

  int sizeLocal = n_body() * (customValues_.elemBufSize(OPERATION_RESTART, NULL, dummy, dummy, dummy) + 4);
  int sizeGlobal = 0, sizeOne = 0;

  // allocate send buffer and pack element data
  // all local elements are in list

  memory->create(sendbuf, sizeLocal, "MultiNodeMeshParallel::writeRestart:sendbuf");
  sizeLocal = 0;
  for (int i = 0; i < n_body(); i++)
  {
    x_bound(xbnd, i);
    sizeOne = customValues_.pushElemToBuffer(i, &(sendbuf[sizeLocal + 4]), OPERATION_RESTART, dummy, dummy, dummy);
    sendbuf[sizeLocal] = static_cast<double>(sizeOne + 4);
    sendbuf[sizeLocal + 1] = xbnd[0];
    sendbuf[sizeLocal + 2] = xbnd[1];
    sendbuf[sizeLocal + 3] = xbnd[2];

    sizeLocal += (sizeOne + 4);
  }

  // gather the per-element data

  sizeGlobal = MPI_Gather0_Vector(sendbuf, sizeLocal, recvbuf, world);

  // write data to file
  if (comm->me == 0)
  {

    // size with 1 extra value (nba)
    int size = (1 + sizeGlobal) * sizeof(double);

    // write size
    fwrite(&size, sizeof(int), 1, fp);

    // write extra value
    fwrite(&nba, sizeof(double), 1, fp);

    // write per-element data
    fwrite(recvbuf, sizeof(double), sizeGlobal, fp);
  }

  // clean up

  memory->destroy(sendbuf);

  if (recvbuf)
    delete[] recvbuf;
}

/* ----------------------------------------------------------------------
   restart functionality - read all required data from restart buffer
   executed on all processes
------------------------------------------------------------------------- */

void MultisphereParallel::restart(double *list)
{
  bool dummy = false;
  int m = 0, nrecv_this;

  int nbody_all_old = static_cast<int>(list[m++]);

  nbody_ = nbody_all_ = 0;

  for (int i = 0; i < nbody_all_old; i++)
  {
    nrecv_this = static_cast<int>(list[m]);

    double *x_bnd = &(list[m + 1]);

    if (domain->is_in_subdomain(x_bnd))
    {

      customValues_.addZeroElement();
      customValues_.deleteRestartElement(nbody_, dummy, dummy, dummy);
      customValues_.popElemFromBuffer(&(list[m + 4]), OPERATION_RESTART, dummy, dummy, dummy);

      nbody_++;
    }
    m += nrecv_this;
  }

  // do initialization tasks

  MPI_Sum_Scalar(nbody_, nbody_all_, world);
  generate_map();
  reset_forces(true);
}
