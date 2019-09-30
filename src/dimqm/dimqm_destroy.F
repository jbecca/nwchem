      subroutine dimqm_destroy(rtdb)
!   Called at end of SCF to destroy DIM/QM arrays
!   Arrays are left allocated until end of calculation
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

      stat = stat .and. ga_destroy(g_DIMxyz)
      stat = stat .and. ga_destroy(g_DIMmuold)
      stat = stat .and. ga_destroy(g_DIMmuind)
     
      if (.not. stat) 
     $  call errquit('dimqm_destroy GA destroy failed', 0, GA_ERR)

    

      end subroutine dimqm_destroy