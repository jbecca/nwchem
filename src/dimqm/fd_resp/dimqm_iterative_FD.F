      subroutine dimqm_iterative_FD(rtdb, g_fld, xyz, nDIM)
c-------------------------------------------------------------------------
c     Calls: qmr_complex2
c     Called by: 
c     Purpose: sets up and calls the qmr solver for FD response
c-------------------------------------------------------------------------
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "dimqm_constants.fh"
#include "dimqm.fh"
c-------------------------- Input variables --------------------------------
      integer rtdb
      integer g_fld
      double precision xyz(3, nDIM)
      integer nDIM
c-------------------------- Local Variables --------------------------------
      integer g_dimwork(7)
      integer i, i3, i0, id, nDIMTypes, nIter, ierr
      integer g_dim_diag_cmplx
      external matvecComplex
      double precision dimtol
      logical stat
      logical ldebug
c
c     Initallize
      id = ga_nodeid()
      call dimqm_get_debug(rtdb, ldebug)
      i3 = 3 * nDIM
!c
      if(.not.rtdb_get(rtdb,'dimqm:ndimtypes',mt_int,1,nDIMTypes))
     $  call errquit('dim_iter2_fd: get nDIMTypes failed',1,RTDB_ERR)

      if(.not.ga_create(mt_dcpl,i3,1,'g_dimwork1',0,0,g_dimwork(1)))
     $ call errquit('g_dimwork1: create failed',0,GA_ERR)
      if(.not.ga_create(mt_dcpl,i3,1,'g_dimwork2',0,0,g_dimwork(2)))
     $ call errquit('g_dimwork2: create failed',0,GA_ERR)
      if(.not.ga_create(mt_dcpl,i3,1,'g_dimwork3',0,0,g_dimwork(3)))
     $ call errquit('g_dimwork3: create failed',0,GA_ERR)
      if(.not.ga_create(mt_dcpl,i3,1,'g_dimwork4',0,0,g_dimwork(4)))
     $ call errquit('g_dimwork4: create failed',0,GA_ERR)
      if(.not.ga_create(mt_dcpl,i3,1,'g_dimwork5',0,0,g_dimwork(5)))
     $ call errquit('g_dimwork5: create failed',0,GA_ERR)
      if(.not.ga_create(mt_dcpl,i3,1,'g_dimwork6',0,0,g_dimwork(6)))
     $ call errquit('g_dimwork6: create failed',0,GA_ERR)
      if(.not.ga_create(mt_dcpl,i3,1,'g_dimwork7',0,0,g_dimwork(7)))
     $ call errquit('g_dimwork7: create failed',0,GA_ERR)

      nIter = MAX(i3, 1000)

      call dimqm_prep_cplx(rtdb,g_dim_diag_cmplx, nDIM, nDIMTypes)
      call dimqm_setlfirst(rtdb,.true.)
      if(.not.rtdb_get(rtdb,'dimqm:dimtol',mt_dbl,1,dimtol))
     $ call errquit('dimqm_iter_FD: get dimtol failed',1,RTDB_ERR)
      call qmr_complex2(rtdb,
     $                  i3,
     $                  g_dimwork,
     $                  nIter,
     $                  dimtol,
     $                  ierr,
     $                  g_fld,
     $                  matvecComplex,
     $                  xyz,
     $                  ldebug,
     $                  g_dim_diag_cmplx,
     $                  nDIM,
     $                  nDIMTypes)

      select case (ierr)
      case (0)
         if (id.eq.0.and.ldebug) then
            write(luout,*)'DIM/QM dipoles converged in ',nIter,'iters'
         endif
      case(2)
         call errquit('error in arguments for QMR cmplx',ierr,GA_ERR)
      case(3)
         call errquit('Iter solver did not converge',ierr,GA_ERR)
      case(4)
         call errquit('Nonsingular system in solver',0,GA_ERR)
      case(5)
         call errquit('A-invariant subspace in solver',0,GA_ERR)
      end select
      stat = .true.
      do i = 1, 7
            stat = stat.and.ga_destroy(g_dimwork(i))
      enddo
      if (.not.stat) 
     $ call errquit('dim_iter_FD:destroy gdimwork failed',0,GA_ERR)
      if (.not.ga_destroy(g_dim_diag_cmplx))
     $ call errquit('dim_iter_FD:destroy dim diag failed',0,GA_ERR)

      end subroutine dimqm_iterative_FD
