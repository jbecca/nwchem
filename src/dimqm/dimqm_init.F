      subroutine dimqm_init(rtdb, calctype)
c       Allocate, initialize, and broadcast persistent arrays
        implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "dimqm.fh"
#include "msgids.fh"

      integer rtdb
      integer calctype
      logical status
c      integer nDim

c      call dimqm_getnDIM(rtdb, nDim)
      

      status = .true.
      if (calctype == 0) then   ! ground state
        status = status .and. (ma_alloc_get
     $            (mt_dbl,3*nDIMAtoms,'dimqm:muind',l_muind, k_muind)) 
        status = status .and. (ma_alloc_get
     $            (mt_dbl,3*nDIMAtoms,'dimqm:muold',l_muold, k_muold)) 
      else if (calctype == 1) then
        status = status .and. (ma_alloc_get
     $            (mt_dcpl,3*nDIMAtoms,'dimqm:muind',l_muind, k_muind)) 
        status = status .and. (ma_alloc_get
     $            (mt_dcpl,3*nDIMAtoms,'dimqm:muold',l_muold, k_muold)) 
      else
        call errquit('dimqm_init: Invalid calctype', 1, MA_ERR)
      endif
      if (.not. status)
     $     call errquit('dimqm_init: MA alloc failed', 1, MA_ERR)

      end