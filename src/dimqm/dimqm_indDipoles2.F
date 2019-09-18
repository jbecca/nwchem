      subroutine dimqm_indDipoles2(rtdb, nbf, geom, basis,
     $                              g_dens, nDIM)
c     Determines the induced dipoles for the given density matrix 
c     This version works for just static response
c     Called from: dimqm_rohf.F
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "dimqm_constants.fh"
#include "geom.fh"
c------------------------ Input Variables --------------------------
      integer rtdb                ! RTDB handle
      integer g_dens              ! perturbed density matrix handle
      integer nbf                 ! number of basis functions
      integer geom                ! geometry handle
      integer basis               ! basis set handle
      integer nDIM
c------------------------ Local variables ---------------------------
      double precision dx, dy, dz
      integer i3, n
      integer l_fld, k_fld
      integer l_dimxyz, k_dimxyz
      integer l_muind, k_muind
      integer l_muold, k_muold
      integer l_eqme, k_eqme
      integer g_dens_comp
      integer alo(3), ahi(3)
      integer blo(2), bhi(2)
      double precision dimErr(3)
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      character*(1) direction(3)
      data direction /'x', 'y', 'z'/
      double precision maxErr, calcErr
      external calcErr
      logical firsttime
      integer id
      logical dimqm_noseed
      logical ldebug
      logical dimqm_seeded
      logical xyz_seeded(6), stat
      double precision dimtol, dimtol0, dimtolxyz(6)
      double precision dsum
      external dsum
c------------------------- Global Arrays -------------------------
      integer g_eqme, g_eqmn, g_muind
c-----------------------------------------------------------------
      id = ga_nodeid()
      call dimqm_get_debug(rtdb, ldebug)
      if(id.eq.0 .and. ldebug) write(luout,*)
     $                  "Starting dimqm_indDipoles2"
      i3 = nDIM * 3
c
      if(.not.ma_push_get(mt_dbl,3*nDIM,'dimqm:dimxyz',l_dimxyz,
     $   k_dimxyz)) call errquit('dimqm_main2: malloc k_dimxyz failed',
     $   1,MA_ERR)
      if(.not.ma_alloc_get(mt_dbl,3*nDIM,'dimqm:eqme',l_eqme,
     $   k_eqme)) call errquit('dimqm_main2: malloc k_eqme failed',
     $   1,MA_ERR)
      if(.not.rtdb_get(rtdb,'dimpar:coords',mt_dbl,
     $      3*nDIM,dbl_mb(k_dimxyz))) call errquit
     $      ('dimmain2: xyz rtdb get failed',1,RTDB_ERR)
c      call dimqm_get_noseed(rtdb, dimqm_noseed)
      if(.not.ga_create(mt_dbl,i3,1,'g_eqme',0,0,g_eqme))
     $   call errquit('dimqm_indDip2: create g_eqme failed',0,GA_ERR)
      if(.not.ga_create(mt_dbl,i3,1,'g_muind',0,0,g_muind))
     $   call errquit('dimqm_indDip2: create g_muind failed',0,GA_ERR)
c
      g_dens_comp = ga_create_atom_blocked(geom,basis,
     $                          'density matrix comp')

      dimErr = ZERO
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf
c
      do n = 1, 3 ! Loop over perturbations
        call ga_zero(g_dens_comp)
        alo(1) = n
        ahi(1) = n 
c
c       Copy current perturbation into g_dens_comp
        call nga_copy_patch('N',g_dens, alo, ahi,
     $                          g_dens_comp, blo, bhi)
c
        call ga_sync() 
c

c       Initizalize arrays
        call ga_zero(g_muind)
        call ga_zero(g_eqme)
        firsttime = .false.
c       Set covergence tolerance
        dimtol = 1.0d-4
        if (.not. rtdb_put(rtdb, 'dimqm:dimtol',mt_dbl,1,dimtol))
     $   call errquit('dim_indDip: put dimtol failed',1,RTDB_ERR)
c       Calculate the QM E-field from the QM density
        call dfill(3*nDIM, 0.0d0, dbl_mb(k_eqme),1)
        call dimqm_EqmE2(rtdb, g_dens_comp, geom, basis,
     $         dbl_mb(k_eqme), dbl_mb(k_dimxyz), nDIM, g_eqme)
c
c       Calculate the dipoles induced by the QM density
        call dimqm_iterative2(rtdb, g_eqme, dbl_mb(k_dimxyz), nDIM)
c
        call ga_init_fence()
        call ga_get(g_eqme, 1, i3, 1, 1, dbl_mb(k_muind), 1)
        call ga_fence()
c       Write induced dipoles to RTDB
        dx = dsum(nDIM, dbl_mb(k_muind), 3)
        dy = dsum(nDIM, dbl_mb(k_muind+1), 3)
        dz = dsum(nDIM, dbl_mb(k_muind+2), 3)
        if(id.eq.0.and.ldebug) then
          write(luout,*) "Total induced dipole moment for "//
     $                  direction(n)//" perturbation"
          write(luout,*) dx, dy, dz
          write(luout,*) ''
        end if
        if(.not.rtdb_put(rtdb, 'dimqm:muind_'//direction(n), mt_dbl,
     $                                 i3, dbl_mb(k_muind)))
     $    call errquit('put dimqm:muind_p failed',1,RTDB_ERR)
      end do ! n = 1, 3
c
c     Clean up arrays
      stat = .true.
      stat = stat .and. ma_free_heap(l_eqme)
      stat = stat .and. ma_chop_stack(l_dimxyz)
      stat = stat .and. ga_destroy(g_eqme)
      stat = stat .and. ga_destroy(g_muind)
      stat = stat .and. ga_destroy(g_dens_comp)
      if(.not.stat)
     $  call errquit('dimqm_indDip2 failed clean up',1,GA_ERR)
c
      end subroutine dimqm_indDipoles2
