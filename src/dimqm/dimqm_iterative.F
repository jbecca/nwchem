      subroutine dimqm_iterative(rtdb, fld, muold, xyz, icmplx, iscmplx,
     $                            dir, nDIM)
c      use constants
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
#include "dimqm_constants.fh"
c
c     Input Variables
      integer rtdb
      integer dtype 
      double precision xyz(3,nDIM)
      double precision fld(3*nDIM, icmplx)
c      double precision muold(nOrder, icmplx)
      double precision muold(3, nDIM, icmplx)
      integer icmplx
      logical iscmplx
      character dir
      logical dimqm_noseed, dimqm_seeded
      integer nDIM, nIter
c
      
c     GA handles
      integer g_work(7)
c
c     Common variables
c      double precision dimtol
c      logical ldebug
c      integer nDIM, nOrder
c
c     Local variables
      integer ierr
      external matvecReal
      external matvecComplex
      double precision time
      logical stat
      integer i, iOrder
      integer id
      double precision err
      integer k_cmplxFld, l_cmplxFld
      integer k_cmplxMuold, l_cmplxMuold
      double precision seed(3*nDIM, icmplx)
      integer g_dim_diag_real
      integer g_dim_diag_cplx
      logical ldebug
      integer nOrder, nDIMTypes
      double precision dimtol

c     Get node ID
      id = ga_nodeid()
      call dimqm_get_debug(rtdb, ldebug)
      if(.not. rtdb_get(rtdb, 'dimpar:ndimtypes', mt_int, 1, nDIMTypes))
     $   call errquit('dimqm_iter: get nDIMTypes failed',1,RTDB_ERR)
      if(id.eq.0.and.ldebug) write(LuOut,*)
     $   "Begin dimqm_iterative"
c
c     Get data type of the external field
      dtype = mt_dbl
      if(iscmplx) dtype = mt_dcpl
c
c     Allocate the work arrays
      nOrder = 3*nDIM
      if(.not. ga_create(dtype, nOrder, 1, 'g_work1', -1, 1,
     $  g_work(1)))
     $  call errquit('g_work1: ga create failed',0,GA_ERR)
      if(.not. ga_create(dtype, nOrder, 1, 'g_work2', -1, 1,
     $  g_work(2)))
     $  call errquit('g_work2: ga create failed',0,GA_ERR)
      if(.not. ga_create(dtype, nOrder, 1, 'g_work3', -1, 1,
     $  g_work(3)))
     $  call errquit('g_work3: ga create failed',0,GA_ERR)
      if(.not. ga_create(dtype, nOrder, 1, 'g_work4', -1, 1,
     $  g_work(4)))
     $  call errquit('g_work4: ga create failed',0,GA_ERR)
      if(.not. ga_create(dtype, nOrder, 1, 'g_work5', -1, 1,
     $  g_work(5)))
     $  call errquit('g_work5: ga create failed',0,GA_ERR)
      if(.not. ga_create(dtype, nOrder, 1, 'g_work6', -1, 1,
     $  g_work(6)))
     $  call errquit('g_work6: ga create failed',0,GA_ERR)
      if(iscmplx) then
        if(.not. ga_create(dtype, nOrder, 1, 'g_work7', -1, 1,
     $    g_work(7)))
     $    call errquit('g_work7: ga create failed',0,GA_ERR)
      end if
c     
c     Some prep before calling the solver
      nIter = MAX(nOrder, 1000)
      if (.not. rtdb_get(rtdb, 'dimqm:dimtol',mt_dbl,1,dimtol))
     $   call errquit('dimqm_iter:get dimtol failed',1,RTDB_ERR)
      err = dimtol
      call dimqm_get_noseed(rtdb, dimqm_noseed)
c      call dimqm_get_seeded(rtdb, dimqm_seeded)

c
c     Call the solver with the proper matrix vector multiplication
c     routine
      time = util_timer()
c      write(luout,*) "noseed?", dimqm_noseed
      if(iscmplx) then
c       
c           Transform fld vector from real to complex
        if(.not.ma_push_get(mt_dcpl, nOrder,'dimqm:cmplxFld',
     $                      l_cmplxFld, k_cmplxFld))
     $   call errquit('dimqm_f2d malloc k_cmplxArr failed', 1, MA_ERR)
c        if(.not.ma_push_get(mt_dcpl, nOrder,'dimqm:cmplxMuold',
c     $                      l_cmplxMuold, k_cmplxMuold))
c     $   call errquit('dimqm_f2d malloc k_cmplxMuold failed', 1, MA_ERR)
        do iOrder = 1, nOrder
          dcpl_mb(k_cmplxFld+iOrder-1) = DCMPLX(fld(iOrder,1),
     $                                          fld(iOrder,2))
c          dcpl_mb(k_cmplxMuold+iOrder-1) = DCMPLX(muold(iOrder,1),
c     $                                            muold(iOrder,2))
        end do
         call dimqm_prep_cplx(rtdb, g_dim_diag_cplx, nDIM,nDIMTypes)
        if(dimqm_noseed) then
          if(id.eq.0 .and. ldebug) write(luout,*) "DIM/QM QMR"
          call qmr_complex(rtdb, nOrder, g_work, nIter, err, ierr, 
     $                     dcpl_mb(k_cmplxFld), matvecComplex, xyz,
     $                     ldebug, g_dim_diag_cplx, nDIM, nDIMTypes)
        else if(dimqm_seeded) then
          if(id.eq.0.and.ldebug) write(luout,*) "Using previous dipoles"
          call seeded_qmr_complex(rtdb, nOrder, g_work, nIter, err,
     $                            ierr, dcpl_mb(k_cmplxFld),
     $                            matvecComplex, xyz, ldebug, 
     $                            dir, g_dim_diag_cplx, nDIM, nDIMTypes)
        else
          if(id.eq.0.and.ldebug) write(luout,*) "DIM/QM Seeding QMR"
          call qmr_seed_complex(rtdb, nOrder, g_work, nIter, err, ierr, 
     $                    dcpl_mb(k_cmplxFld), matvecComplex, xyz,
     $                    ldebug, dir, g_dim_diag_cplx, nDIM, nDIMTypes)
        end if
c
c       Transform fld vector from complex to real
        do iOrder = 1, nOrder
          fld(iOrder, 1) =  REAL(dcpl_mb(k_cmplxFld+iOrder-1))
          fld(iOrder, 2) = AIMAG(dcpl_mb(k_cmplxFld+iOrder-1))
        end do
      else
         call dimqm_prep_real(rtdb, g_dim_diag_real, nDIM,nDIMTypes)
         call dimqm_setlfirst(rtdb,.true.)
c         write(luout,*)'dimqm_noseed',dimqm_noseed
c         write(luout,*)'dimqm_seeded',dimqm_seeded
        if(dimqm_noseed) then
          if(id.eq.0.and.ldebug) write(luout,*) "DIM/QM QMR"
          call qmr_real(rtdb, nOrder, g_work, nIter, err, ierr,
     $                  fld(:,1), matvecReal, xyz,
     $                  ldebug, g_dim_diag_real, nDIM, nDIMTypes)
c           return
        else if(dimqm_seeded) then 
          if(id.eq.0.and.ldebug) write(luout,*) "Using previous dipoles"
          call seeded_qmr_real(rtdb, nOrder, g_work, nIter, err, ierr,
     $                    fld(:,1), matvecReal, xyz, 
     $                    ldebug, dir, g_dim_diag_real, nDIM, nDIMTypes)
        else
          if(id.eq.0.and.ldebug) write(luout,*) "DIM/QM Seeding QMR"
          call qmr_seed_real(rtdb, nOrder, g_work, nIter, err, ierr,
     $               fld(:,1), matvecReal, 
     $               xyz, ldebug, dir, g_dim_diag_real, nDIM, nDIMTypes)
        end if
      end if

      time = util_timer() - time
      select case (ierr)
      case (0)
c       Success
        if(id .eq. 0 .and. ldebug) then
          write(LuOut,*)"DIM/QM dipoles converged in",nIter,"iterations"
          write(LuOut,*) "Took", time, "seconds"
          write(LuOut,*) ''
        end if
      case (2)
        call errquit('Error in arguments for QMR routine',ierr,GA_ERR)
      case (3)
        call errquit('Iterative solver did not converge',ierr,GA_ERR)
      case (4)
        call errquit('Nonsingular system in iterative solver',0,GA_ERR)
      case (5)
        call errquit('A-invarient subspace found in iterative solver',0,
     $                  GA_ERR)
      end select
c
c     Clean up
      do i=1,6 
        stat = ga_destroy(g_work(i))
      end do
      if(iscmplx) then
        stat = ga_destroy(g_work(7))
      end if
      if(iscmplx) then
        if(.not.ma_chop_stack(l_cmplxFld))
     $      call errquit('pop cmplxFld in iterative', 1, MA_ERR)
      end if
      if(iscmplx) then
         if (.not.ga_destroy(g_dim_diag_cplx))
     $      call errquit('dimqm_iter: g_dim_diag_cplx dest?',0,GA_ERR)
      else
         if (.not.ga_destroy(g_dim_diag_real))
     $      call errquit('dimqm_iter: g_dim_diag_real dest?',0,GA_ERR)
      endif

      if(id.eq.0.and.ldebug) write(LuOut,*)
     $   "End dimqm_iterative"
      end subroutine dimqm_iterative
