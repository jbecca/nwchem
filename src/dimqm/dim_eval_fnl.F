c
c
      Subroutine dim_eval_fnl(ptnl,nq,qwght,xyz,qxyz,npert,
     $                        ipm, imag, muind, rtdb)
c
      implicit none
c
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "dimqm_constants.fh"
#include "errquit.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "dimqm.fh"
c
c   Input Variables
      double precision ptnl(nq, npert) ! DIM potential at each quadrature point for each perturbation
      integer nq                       ! Number of quadrature points
      double precision qwght(nq)       ! Quadrature weights
      double precision xyz(3, DIM_nAtoms)    ! DIM atom coordinates
      double precision qxyz(3, nq)     ! Quadrature point coodinates
      integer npert
      integer ipm                      
      integer imag
      double precision muind(3, DIM_nAtoms, npert)  ! Array that holds induced dipoles.  Structured as (component, atom, perturbation)
      integer rtdb
c
c   Local Variables
      integer jj, n, i, idir, id
      double precision r(3), dist, dist3
      double precision screen_mu, screen
      double precision util_erf, scrnFactor
      integer scrnType
c      character*(1) direction(3)
c      character*(1) dpm(2)
c      character*(1) dri(2)
c      data direction /'x', 'y', 'z'/
c      data dpm /'+', '-'/
c      data dri /'r', 'i'/
c
      id = ga_nodeid()
      if (DIM_ldebug.and.(id.eq.0)) then
            write(luout,*)'starting dim_eval_fnl'
            write(luout,*)'npert :', npert
      endif
c
      if (.not. rtdb_get(rtdb, 'dimqm:scrnType', mt_int, 1, scrnType))
     $   call errquit('dim_eval_fnl: scrnType get failed',1,RTDB_ERR)

c   Zero potential array
      call dfill(nq*npert, 0.d0, ptnl, 1)
c      muind = ZERO
      do idir = 1, npert ! Loop over perturbations
c
c        if(id .eq.0) then
c        if(npert > 1) then ! only 1 perturbation for ground state, 3 for response
c          if(ipm .gt. 0) then ! FD Response (real and imaginary, +/-)
c            if(.not.rtdb_get(dimqm_rtdb,'dimqm:muind_'//direction(idir)
c     $                        //'_'//dri(imag)//dpm(ipm), mt_dbl,
c     $                                                  3*nDIM, muind))
c     $      call errquit('get perturbed +/- dipoles failed',1,RTDB_ERR)
c          else ! Static response
c            if(.not.rtdb_get(dimqm_rtdb,'dimqm:muind_'//direction(idir),
c     $                                        mt_dbl,3*nDIM, muind))
c     $      call errquit('get perturbed xyz dipoles failed',1,RTDB_ERR)
c          end if
c        else ! Ground state
c          if(.not.rtdb_get(dimqm_rtdb,'dimqm:muind',
c     $                                        mt_dbl,3*nDIM, muind))
c     $        call errquit('get perturbed dipoles failed',1,RTDB_ERR)
c        endif
c        end if
c       
c        write(luout,*) "muind", id, muind
        if (scrnType .ne. NOSCRN) then
        if (.not. rtdb_get(rtdb,'dimqm:scrnFactor',mt_dbl,1,scrnFactor))
     $   call errquit('dim_eval_fnl:get scrnFactor failed',1,RTDB_ERR)
        endif
        do n = 1, nq ! Loop over quadrature points
          do i = 1, DIM_nAtoms! Loop over DIM atoms
c           Distance between two points
            r(:) = xyz(:,i) - qxyz(:,n) 
            dist = SQRT((r(1)**2 + r(2)**2 + r(3)**2))
c           Screening
            screen_mu = ONE
            select case (scrnType)
            case(ERFSCRN)
              screen = util_erf(scrnFactor * dist)
              screen_mu = screen * screen
            case(EXPSCRN)
              screen = 1.0d0 - EXP(-dist*scrnFactor)
              screen_mu = screen * screen * screen
            case default
              if (dist < 1.0d-12) dist = 1.0d-12
            end select
            dist3 = dist * dist * dist
c           Calculate potential
            ptnl(n,idir) = ptnl(n,idir)
     $                           - screen_mu*muind(1,i,idir)*r(1)/dist3
            ptnl(n,idir) = ptnl(n,idir)
     $                           - screen_mu*muind(2,i,idir)*r(2)/dist3
            ptnl(n,idir) = ptnl(n,idir)
     $                           - screen_mu*muind(3,i,idir)*r(3)/dist3
          end do ! End loop over DIM atoms
        end do ! End loop over quadratrue points
c       Combine with quadrature weights
        do jj = 1, nq
          ptnl(jj,idir) = -ptnl(jj,idir)*qwght(jj)
        enddo
      end do ! End loop over perturbations
      return
      end subroutine dim_eval_fnl
