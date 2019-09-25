      subroutine dimqm_tidy(rtdb)
c     clean up the allocated arrays and variables.
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "dimqm.fh"

      integer rtdb
      logical status


      status = .true.
        status = status .and. (ma_free_heap(l_muind))
        status = status .and. (ma_free_heap(l_muold))
      if (.not. status)
     $     call errquit('dimqm_init: MA free failed', 1, MA_ERR)

      end
