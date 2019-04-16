      subroutine dimqm_f2d2(rtdb, g_etot, g_muind, dimxyz, nDIM)
c-----------------------------------------------------------------------
c  Author: Jeff Becca
c  Calls: dimqm_iterative2
c  Called from: dimqm_main2
c  Purpose: wrapper routine that prepares the fld to be used in solver
c-----------------------------------------------------------------------
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"

      integer rtdb, g_etot, g_muind, nDIM, nDIMTypes, id
      double precision dimxyz(3, nDIM)
      id = ga_nodeid()

c      call dimqm_direct(rtdb, dimxyz, g_etot, nDIM, nDIMTypes)
c  direct solver does not work as of yet. Easy solution was to use 
c  ga_solve, but that ga function requires scalapack. I have not been
c  successful in linking scalapack, so back to iterative it is. 
      
      call dimqm_iterative2(rtdb, g_etot, dimxyz, nDIM)
      
      if(id.eq.0) write(luout,*)'exiting dimqm_f2d2'
      end subroutine dimqm_f2d2
