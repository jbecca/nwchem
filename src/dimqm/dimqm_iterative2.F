      subroutine dimqm_iterative2(rtdb, g_fld, xyz, nDIM)
c-----------------------------------------------------------------------
c  Author: Jeff Becca
c  Calls: qmr_real
c  Called by: dimqm_f2d
c  Purpose: Sets up and calls the qmr solver for the ground state 
c           calculation to solve for DIM dipoles from the given field
c-----------------------------------------------------------------------
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
c---------------------------- Input Variables --------------------------
      integer rtdb
      integer g_fld, g_muind
      double precision xyz(3, nDIM)
      integer nDIM
c---------------------------- Local Variables --------------------------
      integer g_dimwork(6)
      integer id, nDIMTypes, i3, nIter, ierr, i
      logical ldebug, stat
      integer g_dim_diag_real
      external matvecReal
      double precision dimtol

      id = ga_nodeid()
      ldebug = .true.
      if (ldebug.and.id.eq.0) write(luout,*)'starting dimqm_iterative2'
      if(.not.rtdb_get(rtdb, 'dimpar:ndimtypes', mt_int,1,nDIMTypes))
     $   call errquit('dim_iter2: get nDIMTypes failed',1,RTDB_ERR)
      if (ldebug.and.id.eq.0) write(luout,*)'got ndimtypes'
c  Allocate the work GAs
      i3 = 3*nDIM
      if(.not.ga_create(mt_dbl,i3,1,'g_dimwork1',0,0,g_dimwork(1)))
     $   call errquit('g_dimwork1: create failed',0,GA_ERR)
c      call errquit('stop here in iter2')
      if(.not.ga_create(mt_dbl,i3,1,'g_dimwork2',0,0,g_dimwork(2)))
     $   call errquit('g_dimwork2: create failed',0,GA_ERR)
      if(.not.ga_create(mt_dbl,i3,1,'g_dimwork3',0,0,g_dimwork(3)))
     $   call errquit('g_dimwork3: create failed',0,GA_ERR)
      if(.not.ga_create(mt_dbl,i3,1,'g_dimwork4',0,0,g_dimwork(4)))
     $   call errquit('g_dimwork4: create failed',0,GA_ERR)
      if(.not.ga_create(mt_dbl,i3,1,'g_dimwork5',0,0,g_dimwork(5)))
     $   call errquit('g_dimwork5: create failed',0,GA_ERR)
      if(.not.ga_create(mt_dbl,i3,1,'g_dimwork6',0,0,g_dimwork(6)))
     $   call errquit('g_dimwork6: create failed',0,GA_ERR)
      if (ldebug.and.id.eq.0) write(luout,*)'allocated gdimwork arrays'

      nIter = MAX(i3, 1000)

      if (ldebug.and.id.eq.0) write(luout,*)'calling prep real'
      call dimqm_prep_real(rtdb, g_dim_diag_real, nDIM, nDIMTypes)
      if (ldebug.and.id.eq.0) write(luout,*)'back from prep real'
      call dimqm_setlfirst(rtdb,.true.)
      if (ldebug.and.id.eq.0) write(luout,*)'set lfirst'
      if(.not.rtdb_get(rtdb,'dimqm:dimtol',mt_dbl,1,dimtol))
     $   call errquit('dimqm_iter2: get dimtol failed',1,RTDB_ERR)
      if (ldebug.and.id.eq.0) write(luout,*)'calling qmr real'
      call qmr_real2(rtdb, i3, g_dimwork, nIter, dimtol, ierr, g_fld,
     $      matvecReal, xyz, ldebug, g_dim_diag_real, nDIM, nDIMTypes)
      if (ldebug.and.id.eq.0) write(luout,*)'back from qmr real'

      select case (ierr)
      case (0)
c  converged
         if (id.eq.0.and.ldebug) then
            write(luout,*)'DIM/QM dipoles converged in',nIter,'iters'
         endif
      case (2)
         call errquit('error in arguments for QMR real',ierr,GA_ERR)
      case (3)
         call errquit('Iter solver did not converge',ierr,GA_ERR)
      case (4)
         call errquit('Nonsingular system in solver',0,GA_ERR)
      case (5)
         call errquit('A-invariant subspace in solver',0,GA_ERR)
      end select
      stat = .true.
      do i=1,6
         stat = stat.and.ga_destroy(g_dimwork(i))
      enddo
      if (.not.stat) call errquit('dim_iter2:destroy gdimwork failed',
     $   0,GA_ERR)
      if (.not.ga_destroy(g_dim_diag_real)) 
     $   call errquit('dim_iter2:destroy dim diag failed',0,GA_ERR)
      if (ldebug.and.id.eq.0) write(luout,*)'end of dimqm_iterative2'
      end subroutine dimqm_iterative2
