      subroutine dimqm_init(rtdb, calctype)
c       Make sure information in common blocks is same on all processors

        implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "dimqm.fh"

      integer rtdb
      integer calctype
      logical stat
c      integer nDim


      stat = rtdb_parallel(.true.)

      stat = stat .and.
     $       rtdb_get(rtdb,'dimqm:norder',mt_int,1,DIM_nOrder)
      stat = stat .and. 
     $       rtdb_get(rtdb,'dimqm:ldebug',mt_log,1,DIM_ldebug)
      stat = stat .and. 
     $       rtdb_get(rtdb,'dimqm:ndim',mt_int,1,DIM_nAtoms)
      if (.not.stat)
     $      call errquit('DIMQM init rtdb_get 1 failed',1,RTDB_ERR)

      call dimqm_create(rtdb)

      stat = stat .and.
     $       rtdb_get(rtdb,'dimqm:g_DIMxyz',mt_int,1,g_DIMxyz)
      stat = stat .and. 
     $       rtdb_get(rtdb,'dimqm:g_DIMmuind',mt_int,1,g_DIMmuind)
      stat = stat .and.
     $       rtdb_get(rtdb,'dimqm:g_DIMmuold',mt_int,1,g_DIMmuold)
      if (.not.stat)
     $      call errquit('DIMQM init rtdb_get 2 failed',1,RTDB_ERR)
      end