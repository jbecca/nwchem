      subroutine dimqm_fldpot2(rtdb,d_xyz, q_xyz, q_chrg, 
     $                  g_fld, nNuc, nDIM)
      implicit none
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "dimqm_constants.fh"
c  Input Variables
      integer rtdb, g_fld, nNuc, nDIM
      double precision d_xyz(3, nDIM)
      double precision q_xyz(3, nNuc)
      double precision q_chrg(nNuc)
c  Local Variables
      logical ldebug
      integer i, j
      double precision r(3), dist, dist3
      double precision fld(3)
      double precision screen_fld, screen, scrnFactor
      double precision tmpfld(3)
      integer scrnType
      integer nprocs, id
c  START
      ldebug = .true.
      nprocs = ga_nnodes()
      id = ga_nodeid()
      call ga_init_fence() 
      do i = 1, nDIM
         do j = 1, nNuc
c  Get distance between QM atom and DIM atom
         r(:) = d_xyz(:,i) - q_xyz(:,j)
         dist = SQRT((r(1)**2 + r(2)**2 + r(3)**2))
         dist3= dist**3 
         if (.not.rtdb_get(rtdb,'dimqm:scrnType',mt_int,1,scrnType))
     $     call errquit('dimqm_fldpot2: scrnType get failed',1,RTDB_ERR)
         if (scrnType.ne.NOSCRN) then
         if (.not.rtdb_get(rtdb,'dimqm:scrnFactor',mt_dbl,1,scrnFactor))
     $     call errquit('dimqm_fldpot:scrnFactor get failed',1,RTDB_ERR)
         endif

         screen_fld = ONE
         if(scrnType == ERFSCRN) then
           screen     = erf(scrnFactor*dist)
           screen_fld = screen * screen
         else if(scrnType == EXPSCRN) then
           screen = ONE - EXP(-scrnFactor*dist)
           screen_fld = screen * screen * screen
         else
           if(dist < 1.0d-12) dist = 1.0d-12
         end if

c  TODO: fld array not needed and could get large. just do it in
c        tmpfld where it is only double precision (3)
         fld(1) = screen_fld * q_chrg(j) * r(1) / dist3
         fld(2) = screen_fld * q_chrg(j) * r(2) / dist3
         fld(3) = screen_fld * q_chrg(j) * r(3) / dist3
         call ga_acc(g_fld, 3*i-2, 3*i, 1, 1, fld, 1, 1.0d0/nprocs) 
         enddo
      enddo
      call ga_fence()
      end subroutine dimqm_fldpot2
