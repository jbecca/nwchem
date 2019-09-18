      subroutine dimqm_EqmE2(rtdb, g_dens, geom, basis, 
     $               fld, xyz, nDIM, g_eqme)
c-----------------------------------------------------------------------
c  Author: Jeff Becca
c  Called from:   dimqm_main2
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
c--------------------------- Input Variables ---------------------------
      integer rtdb
      integer g_dens
      integer geom
      integer basis
      integer g_eqme
      double precision fld(3, nDIM)
      double precision xyz(3, nDIM) 
      integer nDIM, nOrder
c--------------------------- Local Variables ---------------------------
      logical status
      integer id
      logical ldebug

c      status = rtdb_parallel(.true.)
      id = ga_nodeid()
c      call dimqm_get_debug(rtdb, ldebug)
c  Currently we call hnd_elfcon for getting the e- integral contribution
c  to the DIM fields, which uses malloc array, but we want the global 
c  arrays for the dipoles and such for large DIM systems, so we copy
c  into global array after solving. 
      ldebug = .true.
c      if (ldebug.and.id.eq.0) write(luout,*)'Starting EqmE2'
      nOrder = 3*nDIM
      call hnd_elfcon_0(basis, geom, g_dens, xyz, nDIM, 
     $   fld, 1, .false., .false.)
      call ga_put(g_eqme, 1, nOrder, 1, 1, fld, 1)
      call ga_sync()
c      call ga_print(g_eqme)
c      if (ldebug.and.id.eq.0) write(luout,*)'Exiting EqmE2'
      end subroutine dimqm_EqmE2


