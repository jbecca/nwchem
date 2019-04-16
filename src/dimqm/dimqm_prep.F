c      subroutine dimqm_prep(rtdb, om, nDIM)
cc
cc     Author: Jeff Becca
cc     Called from: dimqm_main.F cphf_solve2.F cphf_solve3.F 
cc                  aoresponse_driver.F aoresponse_driver_new.F
cc                  dimqm_polarizability.F
cc                  
cc
cc     Subroutines called: None
cc
cc     Builds the interaction matrix for the DIM system.   Called
cc     at the start of every DIM/QM calculation.
cc
cc     Currently stripped down for just an iterative solver. If a 
cc     direct solver or CMM is implemented, prep will need to be altered
cc
c      implicit none
c#include "errquit.fh"
c#include "inp.fh"
c#include "rtdb.fh"
c#include "stdio.fh"
c#include "nwc_const.fh"
c#include "mafdecls.fh"
c#include "global.fh"
c#include "dimqm_constants.fh"
c#include "dimqm.fh"
cc
cc     Input variables
cc
c      integer rtdb        ! RTDB handle
c      double precision om ! Frequency, zero or negative means static
cc
cc     Local variables
c      integer nDIM
c      integer stat
c      integer nLength
c      integer g_matrix
c      integer g_indx
c      integer m, i, itype
c      integer mLo, mUp, mOnNode
c      double precision polr
c      double complex   polc
c      double complex   die
c      integer tx(nDIM) ! Type index
c      double precision sPol(nDIMTypes)
c      double precision bound(nDIMTypes)
c      double precision drudeterms(2, nDIMTypes)
c      double complex Drude
c      double complex ExperimentalDielectric
c      integer id
c      double complex expdies(nDIMTypes)
c      logical expdiesRead(nDIMTypes)
c      logical hasExpDie(nDimTypes)
c      integer g_dim_diag
c      logical lfirst
c      logical ldebug
cc
cc     Common variables used from dimqm.fh
cc
cc      integer alg, g_DIM_diag
cc      logical ldebug, lpim
cc      integer nOrder, nDIM, nDIMTypes
cc
c      expdies = ZERO_C
c      expdiesRead = .false.
c      id = ga_nodeid()
c      if(id.eq.0 .and. ldebug) write(luout,*) "Start dimqm_prep"
c      call pphilo(id, nDIM, mLo, mUp, mOnNode)
cc
cc   Pull needed info from RTDB
c      if(.not. rtdb_get(rtdb,'dimpar:sPol', mt_dbl, nDIMTypes, sPol))
c     $  call errquit('dimpar:sPol rtdb get failed', 0, RTDB_ERR)
cc
c      if(.not. rtdb_get(rtdb,'dimpar:typeindex', mt_int, nDIM, tx))
c     $  call errquit('dimpar:type index rtdb get failed', 0, RTDB_ERR)
cc
c      if(.not. rtdb_get(rtdb,'dimpar:bound', mt_dbl, nDIMTypes, bound))
c     $  call errquit('dimpar:bound rtdb get failed', 0, RTDB_ERR)
cc
c      if(.not. rtdb_get(rtdb,'dimpar:drude', mt_dbl, nDIMTypes*2,
c     $                                    drudeterms))
c     $  call errquit('dimpar:drude rtdb get failed', 0, RTDB_ERR)
c      if(.not. rtdb_get(rtdb,'dimpar:hasexp', mt_log, nDIMTypes,
c     $                                    hasExpDie))
c     $  call errquit('dimpar:hasexp rtdb get failed', 0, RTDB_ERR)
c
cc      if (ldebug) write(luout,*) 'values grabbed from rtdb'
cc      if (ldebug) write(luout,*) 'spol',sPol
cc      if (ldebug) write(luout,*) 'tx',tx
cc      if (ldebug) write(luout,*) 'bound',bound
cc      if (ldebug) write(luout,*) 'drude',drudeterms
cc      if (ldebug) write(luout,*) 'hasexpdie',hasExpDie
c
c     
cc
cc   =========================================
cc    Allocate GAs for the interaction matrix
cc   =========================================
cc
cC      if (ldebug) write(luout,*)'alg',alg
cC      if (ldebug) write(luout,*)'om',om
c
ccc     Exact Solver
cc      if(alg ==  0) then   
cc        nLength = nOrder*(nOrder+1)/2
ccc
ccc        Complex matrix allocation
cc        if(om > ZERO) then   
cc          if(.not.ga_create(mt_dcpl, nLength, 1, 'dimMatrix_FD',
cc     $                      -1, 1, g_matrix))
cc     $      call errquit('dimMatrix_FD: ga_create failed',9,GA_ERR)
ccc
ccc        DP matrix allocation
cc        else
cc          if(.not.ga_create(mt_dbl, nLength, 1, 'dimMatrix_S', -1, 1, 
cc     $                      g_matrix))
cc     $      call errquit('dimMatrix_S: ga_create failed',9,GA_ERR)
cc        end if
cc      if (ldebug) write(luout,*)'dimmatrix allocated'
ccc
ccc        Allocate index array
cc        if(.not. ga_create(mt_int, nOrder, 1, 'dimMatrix_indx', -1, 1,
cc     $                     g_indx))
cc     $    call errquit('dimMatrix_indx: ga_create failed',0,GA_ERR)
cc        if (ldebug) write(luout,*) 'dim index matrix allocated'
ccc
ccc        Zero both arrays
cc        call ga_zero(g_matrix)
cc        call ga_zero(g_indx)
cc        if (ldebug) write(luout,*) 'arrays zerod'
ccc
ccc        Place handles into RTDB
cc        if(.not. rtdb_put(rtdb,'dimqm:g_matrix',mt_int,1,g_matrix))
cc     $    call errquit('g_matrix: rtdb put failed',0,RTDB_ERR)
cc        if(.not. rtdb_put(rtdb,'dimqm:g_indx',mt_int,1,g_indx))
cc     $    call errquit('g_indx: rtdb put failed',0,RTDB_ERR)
cc        if (ldebug) write(luout,*) 'arrays handles placed in rtdb'
ccc
cc      Iterative Solver
cc      else
cc        Complex diagonal allocation
c      if(om > ZERO) then
c         if (ldebug) write(luout,*) 'dim diagFD about to be allocated'
c         if (ldebug) write(luout,*) 'norder',norder
c         if(.not. ga_create(mt_dcpl, nOrder, 1, 'diagFD',
c     $                       -1, -1, g_DIM_diag))
c     $      call errquit('dimqm:alphainvFD: create failed',0,GA_ERR)
cc      if (om > zero) call errquit('i got here in prep',0,GA_ERR)
c         if (ldebug) write(luout,*) 'dim diagFD allocated'
c         call ga_zero(g_DIM_diag)
cc        DP diagonal allocation
c      else
c         if(.not. ga_create(mt_dbl, nOrder, 1, 'dimqm:diagS',
c     $                       -1, -1, g_DIM_diag))
c     $    call errquit('dimqm:alphainvFD: create failed',9,GA_ERR)
c         call ga_zero(g_DIM_diag)
c      end if
c      if (ldebug) write(luout,*) 'complex diag arrays allocated'
cc
cc        Zero the diagonal
c        if (ldebug) write(luout,*) 'complex diag arrays zerod'
cc      end if
cc
cc   =========================================================================
cc   Build the interaction matrix (only the diagonal for the iterative solver) 
cc   =========================================================================
cc
cc      Iterative solver
c      if(alg > 0) then
c        call ga_init_fence()
c        if (ldebug) write(luout,*) 'ga_init_fence called'
cc
cc        Frequency dependent calculation
c        if(om > ZERO) then
cc          Loop over DIM atoms
c          do m = mLo, mUp 
c            itype = tx(m)
c            die = ExperimentalDielectric(rtdb, itype, om, hasExpDie,
c     $                                   expDiesRead, expDies)
c            if (die /= ZERO_C) then
c              ! Size correction
c            else
c              die = CMPLX(bound(itype), ZERO)
c              if(drudeterms(1,itype) > ZERO) then ! Subtract Drude term if given
c              die = die -
c     $              Drude(om, drudeterms(1,itype), drudeterms(2,itype))
c              endif
c            end if
c            polc = ONE_C/(((die - ONE_C) / (die + TWO_C)) * sPol(itype))
c            i = 3*(m-1)
c            call ga_put(g_DIM_diag, i+1, i+1, 1, 1, polc, 1)
c            call ga_put(g_DIM_diag, i+2, i+2, 1, 1, polc, 1)
c            call ga_put(g_DIM_diag, i+3, i+3, 1, 1, polc, 1)
c          end do
c          if (ldebug) write(luout,*) 'ga_put polc in dim diag finished'
cc
cc        Static calculation
c        else
cc          Loop over DIM atoms
c          do m = mLo, mUp
c            itype = tx(m)
c            i = 3*(m-1)
c            polr = ONE / sPol(itype)
c            call ga_put(g_DIM_diag, i+1, i+1, 1, 1, polr, 1)
c            call ga_put(g_DIM_diag, i+2, i+2, 1, 1, polr, 1)
c            call ga_put(g_DIM_diag, i+3, i+3, 1, 1, polr, 1)
c          end do
c        end if
c        call ga_fence()
cc
cc      Exact Solver
c      else
c        call errquit("dimqm direct solver not implemented", 1,
c     $               INPUT_ERR)
c      end if
cc      dimqm_on = .true.
c      lfirst   = .true.
c
c      call dimqm_setlfirst(rtdb, lfirst)
c      
c
c      if (ldebug) write(luout,*)'exiting dimqm_prep'
c      end subroutine dimqm_prep

      function Drude (omega, plasma, damp) result (drd)
      implicit none
#include "dimqm_constants.fh"
c
c   =====================================
c   Calculates the complex Drude function
c   =====================================
c
      double precision omega
      double precision plasma
      double precision damp
      double complex drd

      drd = plasma**2 / (omega * (omega + I_C * damp))
      end function Drude

      function ExperimentalDielectric(rtdb, atomtype, om, hasExpDie,
     $                      expDieRead, expDies, nDIMTypes) result (die)
      implicit none
#include "rtdb.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "dimqm_constants.fh"
#include "errquit.fh"

      integer rtdb
      integer atomtype
      double precision om
      double complex die
      logical hasExpDie(nDIMTypes)
      double complex expDies(nDIMTypes)
      logical expDieRead(nDIMTypes)
      
      integer iexp, nDIMTypes
      double precision rawexpdie(3, 10000)
      double precision splined(2, 10000)
      double complex interpolate_complex1

      if (.not. hasExpDie(atomtype)) then
        die = ZERO_C
      else if (expDieRead(atomtype)) then
        die = expDies(atomtype)
      else
        if(.not.rtdb_get(rtdb,'dimpar:iexp'//char(atomtype),
     $                   mt_int, 1, iexp))
     $      call errquit('iexp: rtdb get failed', atomtype, RTDB_ERR)
        if(.not.rtdb_get(rtdb,'dimpar:expdata'//char(atomtype), mt_dbl,
     $                   3*iexp, rawexpdie(1:3,1:iexp)))
     $      call errquit('expdata: rtdb get failed', atomtype, RTDB_ERR)
        if(.not.rtdb_get(rtdb,'dimpar:splinedata'//char(atomtype),
     $                   mt_dbl, 2*iexp, splined(1:2,1:iexp)))
     $     call errquit('splinedata: rtdb get failed',atomtype,RTDB_ERR)
        expDies(atomtype) = interpolate_complex1(iexp,
     $                                           rawexpdie(1,1:iexp),
     $                                           rawexpdie(2,1:iexp),
     $                                           rawexpdie(3,1:iexp),
     $                                           splined(1,1:iexp),
     $                                           splined(2,1:iexp),
     $                                           HART2NM / om)
        expDieRead(atomtype) = .true.
        die = expDies(atomtype)
      end if
      end function ExperimentalDielectric


      subroutine dimqm_prep_real(rtdb, g_dim_diag_real, nDIM, nDIMTypes)
c
c     Author: Jeff Becca
c
c     Called from: dimqm_iterative
c
c     Subroutines called: None
c
c     Builds the interaction matrix for the DIM system.  
c
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "dimqm_constants.fh"
c
c     Input variables
c
      integer rtdb        ! RTDB handle
      integer g_dim_diag_real    ! real diag of int. mat. 
      integer nDIM, nDIMTypes
c
c     Local variables
      integer stat
      integer nLength
      integer g_matrix
      integer g_indx
      integer m, i, itype
      integer mLo, mUp, mOnNode
      double precision polr
      double complex   polc
      double complex   die
      integer tx(nDIM) ! Type index
      double precision sPol(nDIMTypes)
      double precision bound(nDIMTypes)
      double precision drudeterms(2, nDIMTypes)
      double complex Drude
      double complex ExperimentalDielectric
      integer id
      double complex expdies(nDIMTypes)
      logical expdiesRead(nDIMTypes)
      logical hasExpDie(nDimTypes)
      integer g_dim_diag2
      logical lfirst
      logical ldebug
      integer nOrder
c
c     Common variables used from dimqm.fh
c
c      integer alg, 
c      logical ldebug, lpim
c      integer nOrder, nDIM, nDIMTypes
c
      call dimqm_get_debug(rtdb, ldebug)
      nOrder = 3*nDIM
      expdies = ZERO_C
      expdiesRead = .false.
      id = ga_nodeid()
      if(id.eq.0 .and. ldebug) write(luout,*) "Start dimqm_prep_real"
c      call pphilo(id, nDIM, mLo, mUp, mOnNode)
c
c   Pull needed info from RTDB
      if(.not. rtdb_get(rtdb,'dimpar:sPol', mt_dbl, nDIMTypes, sPol))
     $  call errquit('dimpar:sPol rtdb get failed', 0, RTDB_ERR)
      if(.not. rtdb_get(rtdb,'dimpar:typeindex', mt_int, nDIM, tx))
     $  call errquit('dimpar:type index rtdb get failed', 0, RTDB_ERR)
      if(.not. rtdb_get(rtdb,'dimpar:bound', mt_dbl, nDIMTypes, bound))
     $  call errquit('dimpar:bound rtdb get failed', 0, RTDB_ERR)
      if(.not. rtdb_get(rtdb,'dimpar:drude', mt_dbl, nDIMTypes*2,
     $                                    drudeterms))
     $  call errquit('dimpar:drude rtdb get failed', 0, RTDB_ERR)
      if(.not. rtdb_get(rtdb,'dimpar:hasexp', mt_log, nDIMTypes,
     $                                    hasExpDie))
     $  call errquit('dimpar:hasexp rtdb get failed', 0, RTDB_ERR)

c      if (ldebug) write(luout,*) 'values grabbed from rtdb'
c      if (ldebug) write(luout,*) 'spol',sPol
c      if (ldebug) write(luout,*) 'tx',tx
c      if (ldebug) write(luout,*) 'bound',bound
c      if (ldebug) write(luout,*) 'drude',drudeterms
c      if (ldebug) write(luout,*) 'hasexpdie',hasExpDie
c
c   =========================================
c    Allocate GAs for the interaction matrix
c   =========================================
c

c        DP diagonal allocation
      if(.not. ga_create(mt_dbl, nOrder, 1, 'dimqm:diagS',
     $                       nOrder, -1, g_DIM_diag_real))
     $    call errquit('dimqm:alphainvFD: create failed',9,GA_ERR)

      if (ldebug.and.id.eq.0) write(luout,*) 'diagS arrays allocated'
c
c        Zero the diagonal
      call ga_zero(g_DIM_diag_real)
      if (ldebug.and.id.eq.0) write(luout,*) 'diags arrays zerod'
c
c   =========================================================================
c   Build the interaction matrix (only the diagonal for the iterative solver) 
c   =========================================================================
c
c      Iterative solver
c      if(alg > 0) then
        call ga_init_fence()
        if (ldebug.and.id.eq.0) write(luout,*) 'ga_init_fence called'
c
c          Loop over DIM atoms
          do m = 1, nDIM
            itype = tx(m)
            i = 3*(m-1)
            polr = ONE / sPol(itype)
            call ga_put(g_DIM_diag_real, i+1, i+1, 1, 1, polr, 1)
            call ga_put(g_DIM_diag_real, i+2, i+2, 1, 1, polr, 1)
            call ga_put(g_DIM_diag_real, i+3, i+3, 1, 1, polr, 1)
          end do
        call ga_fence()
c
c      Exact Solver
c      else
c        call errquit("dimqm direct solver not implemented", 1,
c     $               INPUT_ERR)
c      end if
c      dimqm_on = .true.
      lfirst   = .true.

      call dimqm_setlfirst(rtdb, lfirst)

      if (ldebug.and.id.eq.0) write(luout,*)'exiting dimqm_prep_real'
      end subroutine dimqm_prep_real



      subroutine dimqm_prep_cplx(rtdb, g_dim_diag_cplx, nDIM, nDIMTypes)
c
c     Author: Jeff Becca
c
c     Called from: dimqm_iterative
c
c     Subroutines called: None
c
c     Builds the interaction matrix for the DIM system.  
c
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "dimqm_constants.fh"
c
c     Input variables
c
      integer rtdb        ! RTDB handle
c      integer om          ! omega, frequency
      integer g_dim_diag_cplx    ! real diag of int. mat. 
      integer nDIM, nDIMTypes
c
c     Local variables
      integer stat
      integer nLength
      integer g_matrix
      integer g_indx
      integer m, i, itype
      integer mLo, mUp, mOnNode
      double precision polr
      double complex   polc
      double complex   die
      integer tx(nDIM) ! Type index
      double precision sPol(nDIMTypes)
      double precision bound(nDIMTypes)
      double precision drudeterms(2, nDIMTypes)
      double complex Drude
      double complex ExperimentalDielectric
      integer id
      double complex expdies(nDIMTypes)
      logical expdiesRead(nDIMTypes)
      logical hasExpDie(nDimTypes)
      integer g_dim_diag2
      double precision om
      logical lfirst
      logical ldebug
      integer nOrder
c
c     Common variables used from dimqm.fh
c
c      integer alg, 
c      logical ldebug, lpim
c      integer nOrder, nDIM, nDIMTypes
c

      call dimqm_get_debug(rtdb, ldebug)
      nOrder = 3*nDIM
      expdies = ZERO_C
      expdiesRead = .false.
      id = ga_nodeid()
      if(id.eq.0 .and. ldebug) write(luout,*) "Start dimqm_prep_cplx"
c      call pphilo(id, nDIM, mLo, mUp, mOnNode)
c
c   Pull needed info from RTDB
      if(.not. rtdb_get(rtdb,'dimpar:sPol', mt_dbl, nDIMTypes, sPol))
     $  call errquit('dimpar:sPol rtdb get failed', 0, RTDB_ERR)
      if(.not. rtdb_get(rtdb,'dimpar:typeindex', mt_int, nDIM, tx))
     $  call errquit('dimpar:type index rtdb get failed', 0, RTDB_ERR)
      if(.not. rtdb_get(rtdb,'dimpar:bound', mt_dbl, nDIMTypes, bound))
     $  call errquit('dimpar:bound rtdb get failed', 0, RTDB_ERR)
      if(.not. rtdb_get(rtdb,'dimpar:drude', mt_dbl, nDIMTypes*2,
     $                                    drudeterms))
     $  call errquit('dimpar:drude rtdb get failed', 0, RTDB_ERR)
      if(.not. rtdb_get(rtdb,'dimpar:hasexp', mt_log, nDIMTypes,
     $                                    hasExpDie))
     $  call errquit('dimpar:hasexp rtdb get failed', 0, RTDB_ERR)
      if(.not. rtdb_get(rtdb, 'dimqm:omega', mt_dbl, 1, om))
     $  call errquit('dimqm:omega rtdb get failed', 0, RTDB_ERR)
c   =========================================
c    Allocate GAs for the interaction matrix
c   =========================================


c        Complex diag allocation
      if(.not. ga_create(mt_dcpl, nOrder, 1, 'dimqm:diagFD',
     $                       nOrder, -1, g_DIM_diag_cplx))
     $    call errquit('dimqm:alphainvFD: create failed',9,GA_ERR)

      if (ldebug.and.id.eq.0) write(luout,*) 'diagFD arrays allocated'
c
c        Zero the diagonal
      call ga_zero(g_DIM_diag_cplx)
      if (ldebug.and.id.eq.0) write(luout,*) 'diagFD arrays zerod'
c
c   =========================================================================
c   Build the interaction matrix (only the diagonal for the iterative solver) 
c   =========================================================================
c
c      Iterative solver
c      if(alg > 0) then
        call ga_init_fence()
        if (ldebug.and.id.eq.0) write(luout,*) 'ga_init_fence called'
c          Loop over DIM atoms
c          do m = mLo, mUp
          do m = 1, nDIM
            itype = tx(m)
            die = ExperimentalDielectric(rtdb, itype, om, hasExpDie,
     $                                  expDiesRead, expDies, nDIMTypes)
            if (die /= ZERO_C) then
               ! size correction? at some point?
            else
               die = CMPLX(bound(itype), ZERO)
               if(drudeterms(1,itype) > ZERO) then ! subtract drude term given
               die = die - 
     $               Drude(om, drudeterms(1,itype), drudeterms(2,itype))
               endif
            endif
            polc = ONE_C/(((die - ONE_C) / (die + TWO_C)) * sPol(itype))
            i = 3*(m-1)
            call ga_put(g_DIM_diag_cplx, i+1, i+1, 1, 1, polc, 1)
            call ga_put(g_DIM_diag_cplx, i+2, i+2, 1, 1, polc, 1)
            call ga_put(g_DIM_diag_cplx, i+3, i+3, 1, 1, polc, 1)
          end do
        call ga_fence()
c
c      Exact Solver
c      else
c        call errquit("dimqm direct solver not implemented", 1,
c     $               INPUT_ERR)
c      end if
c      dimqm_on = .true.
      lfirst   = .true.

      call dimqm_setlfirst(rtdb, lfirst)

      if (ldebug.and.id.eq.0) write(luout,*)'exiting dimqm_prep_cplx'
      end subroutine dimqm_prep_cplx

