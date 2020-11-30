      subroutine dimqm_eqme_FD(rtdb, g_dens, geom, basis, 
     $               xyz, g_eqme)
c-----------------------------------------------------------------------
c  Author: Jeff Becca
c  Called from:   dimqm_FD_cycle
c  Calls:   hnd_elfcon
c  Purpose: This routine calculates the DIM field contribution from the
c           QM electrons.
c-----------------------------------------------------------------------
      implicit none
#include "rtdb.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
#include "dimqm.fh"
c--------------------------- Input Variables ---------------------------
      integer rtdb
      integer g_dens
      integer geom
      integer basis
      integer g_eqme
      double precision xyz(3, DIM_nAtoms) 
c--------------------------- Local Variables ---------------------------
      integer nOrder
      logical stat
      integer id
      logical ldebug
      double precision fld(3, DIM_nAtoms)
      integer l_fld, k_fld

      id = ga_nodeid()
      ldebug = .true.
      nOrder = 3*DIM_nAtoms
      if(.not.ma_alloc_get(mt_dbl,3*DIM_nAtoms,'dimqm:eqme',l_fld,
     $    k_fld)) call errquit('dimqm: malloc k_fld failed',1,MA_ERR)
      call hnd_elfcon_0(basis, geom, g_dens, xyz, DIM_nAtoms, 
     $   dbl_mb(k_fld), 1, .false., .false.)
      call ga_put(g_eqme, 1, nOrder, 1, 1, dbl_mb(k_fld), 1)
      call ga_sync()
      stat = .true.
      stat = stat.and.ma_free_heap(l_fld)
      end subroutine dimqm_eqme_FD


