      subroutine dimqm_create(rtdb)
*
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "dimqm.fh"

      integer rtdb
      logical stat

      stat = .true.

      stat = stat .and. 
     $  ga_create(mt_dbl,3*DIM_nAtoms,1,'g_DIMxyz', 0, 0, g_DIMxyz)
      stat = stat .and.
     $   rtdb_put(rtdb,'dimqm:g_DIMxyz',mt_int,1,g_DIMxyz)
      stat = stat .and.
     $  ga_create(mt_dbl,3*DIM_nAtoms,1,'g_DIMmuind', 0, 0, g_DIMmuind)
      stat = stat .and.
     $   rtdb_put(rtdb,'dimqm:g_DIMmuind',mt_int,1,g_DIMmuind)
      stat = stat .and.
     $  ga_create(mt_dbl,3*DIM_nAtoms,1,'g_DIMmuold', 0, 0, g_DIMmuold)
      stat = stat .and.
     $   rtdb_put(rtdb,'dimqm:g_DIMmuold',mt_int,1,g_DIMmuold)
     
      if (.not. stat) 
     $  call errquit('dimqm_create GA create failed', 0, GA_ERR)

    

      end subroutine dimqm_create