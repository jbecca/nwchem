      subroutine dimqm_FD_cycle(    rtdb,
     $                              geom,
     $                              basis,
     $                              nbf,
     $                              ncomp, 
     $                              g_dens)
c  ==================================================================
c  Calculates and adds the frequency-dependent DIM potential to the 
c  response Fock matricies (real and imaginary).  Requires knowledge
c  of both the real and imaginary vectors simultaneously which is
c  why this has to be located here, unlike the static routine.
c  ==================================================================
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "dimqm_constants.fh"
#include "dimqm.fh"
#include "geom.fh"
#include "crohf.fh"
c #include "cscf.fh"
c
c     Input Variables
      integer rtdb, geom, basis, nbf
      integer ncomp
      integer g_dens(4)
c
      integer i3, idir
      integer id, n, ipm
      integer g_dens_comp_i, g_dens_comp_r,g_dens_comp_tot
      integer g_fld1, g_fld2, g_fldcmplx
      integer g_fld_tot
      logical stat
      integer dims(2), chunk(2)
      integer alo(3), ahi(3)
      integer blo(2), bhi(2)
      integer mlo, mhi
      integer nlo(2), nhi(2) 
      integer g_mu_r, g_mu_i
      integer g_vdimr, g_vdimi
      integer d(3), dchunk(3)
      integer l_dimxyz, k_dimxyz
c
      id = ga_nodeid()
      i3 = DIM_nAtoms * 3

      dims(1) = nbf
      dims(2) = nbf
      chunk(1) = -1
      chunk(2) = -1

      !NOTE: If this whole thing is written in terms of (ipm) loop, then
      ! can basically copy the dimqm_lclfld2 routine inside of (ipm) loop
      
      if(.not.ga_create(mt_dcpl,i3,3,'g_fld_tot',0,0,g_fld_tot))
     $   call errquit('dim_fd_cycle: create g_fld_tot failed',0,GA_ERR)
      if(.not.nga_create(mt_dbl,2,dims,'dimFD dens_r',
     $   chunk,g_dens_comp_r)) call errquit(
     $   'dimqm_fd: create dens_comp_r failed',0,GA_ERR)
      if(.not.nga_create(mt_dbl,2,dims,'dimFD dens_i',
     $   chunk,g_dens_comp_i)) call errquit(
     $   'dimqm_fd: create dens_comp_i failed',0,GA_ERR)
      if(.not.nga_create(mt_dcpl,2,dims,'dimFD dens_tot',
     $   chunk,g_dens_comp_tot)) call errquit(
     $   'dimqm_fd: create dens_comp_i failed',0,GA_ERR)
      if(.not.ga_create(mt_dbl,DIM_nOrder,1,'g_fld1',-1,-1,g_fld1))
     $   call errquit('dimqm_fd: create g_fld1 failed',0,GA_ERR)
      if(.not.ga_create(mt_dbl,DIM_nOrder,1,'g_fld2',-1,-1,g_fld2))
     $   call errquit('dimqm_fd: create g_fld2 failed',0,GA_ERR)
      if(.not.ga_create(mt_dcpl,DIM_nOrder,1,'g_fldcmplx',
     $      -1,-1,g_fldcmplx)) call errquit(
     $      'dimqm_fd: create g_fldcmplx failed',0,GA_ERR)
      if(.not.ma_push_get(mt_dbl,DIM_nOrder,'dimqm:dimxyz',l_dimxyz,
     $   k_dimxyz)) call errquit('dimqm: malloc k_dimxyz failed',
     $   1, MA_ERR)
      if(.not.rtdb_get(rtdb,'dimqm:coords',mt_dbl,i3,dbl_mb(k_dimxyz)))
     $  call errquit('dim_fd_cycle: rtdb get failed',1,RTDB_ERR)
      if(.not.ga_create(mt_dbl,DIM_nOrder,3,'g_mu_r',
     $      0,0,g_mu_r)) call errquit(
     $      'dimqm_fd: create g_mu_r failed',0,GA_ERR)
      if(.not.ga_create(mt_dbl,DIM_nOrder,3,'g_mu_i',
     $      0,0,g_mu_i)) call errquit(
     $      'dimqm_fd: create g_mu_i failed',0,GA_ERR)
      call ga_zero(g_fld_tot)
      call ga_zero(g_mu_r)
      call ga_zero(g_mu_i)
      d(1)   = 3
      d(2)   = nbf
      d(3)   = nbf
      dchunk(1) = d(1)
      dchunk(2) = -1
      dchunk(3) = -1

      if(.not.nga_create(mt_dbl,3,d,'dimpot_r',dchunk,
     $      g_vdimr)) call errquit(
     $      'dimqm_fd_cycle: create g_dimr+ failed',1,GA_ERR)
      if(.not.nga_create(mt_dbl,3,d,'dimpot_i',dchunk,
     $      g_vdimi)) call errquit(
     $      'dimqm_fd_cycle: create g_dimr+ failed',1,GA_ERR)

      alo(2) = 1 
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf

      mlo = 1
      mhi = DIM_nOrder
      nlo(2) = 1
      nhi(2) = DIM_nOrder

      do ipm = 1, 2 ! loop over perturbation sign
            do n = 1, 3 ! loop over direction
              alo(1) = n
              ahi(1) = n
              nlo(1) = n
              nhi(1) = n
              call ga_zero(g_dens_comp_r)
              call ga_zero(g_dens_comp_i)
              call ga_zero(g_fld1)
              call ga_zero(g_fld2)
              call ga_zero(g_fldcmplx)
              call nga_copy_patch('n', g_dens(ipm), alo, ahi,
     $                               g_dens_comp_r, blo, bhi)
              call nga_copy_patch('n', g_dens(ipm+2), alo, ahi,
     $                                 g_dens_comp_i, blo, bhi)
              call dimqm_EqmE_FD(
     $                              rtdb,
     $                              g_dens_comp_r,
     $                              geom,
     $                              basis,
     $                              dbl_mb(k_dimxyz),
     $                              g_fld1)
              call dimqm_EqmE_FD(
     $                              rtdb,
     $                              g_dens_comp_i,
     $                              geom,
     $                              basis,
     $                              dbl_mb(k_dimxyz),
     $                              g_fld2)
            ! dimqm iterative FD needs to be called here for each 
            ! perturbation direction and plus/ minus 
            ! to make it match how local field routine calls dimqm_iterative_FD
            ! routine convert_d2z is segfaulting on chopping stack. 
            ! trying writing own version of that routine here to 
            ! make our field vector complex
              call dim_convert_d2z(DIM_nOrder,
     $                             g_fld1, 
     $                             g_fld2, 
     $                             g_fldcmplx)
              call dimqm_iterative_FD(
     $                              rtdb,
     $                              g_fldcmplx,
     $                              dbl_mb(k_dimxyz),
     $                              DIM_nAtoms)
              call ga_copy_patch('n', g_fldcmplx, 1, i3, 1, 1,
     $                                g_fld_tot,  1, i3, n, n)
            enddo !n 
            call convert_z2d(g_fld_tot, g_mu_r, g_mu_i)
            call ga_print(g_mu_r)
            call dim_grid_FD_gen(   rtdb,
     $                        g_vdimr,
     $                        g_mu_r,
     $                        1,
     $                        1)
            call dim_grid_FD_gen(   rtdb,
     $                        g_vdimi,
     $                        g_mu_i,
     $                        1,
     $                        1)
            alo(2) = 1
            ahi(2) = nbf
            alo(3) = 1
            ahi(3) = nbf
            call ga_print(g_vdimr)
            call ga_print(g_dens(ipm))
            do idir = 1, 3
                  alo(1) = idir
                  ahi(1) = idir
                  call nga_add_patch(ONE, g_dens(ipm), alo, ahi,
     $                               -ONE, g_vdimr, alo, ahi,
     $                                    g_dens(ipm), alo, ahi)
                  call nga_add_patch(ONE, g_dens(ipm+2), alo, ahi,
     $                               -ONE, g_vdimi, alo, ahi,
     $                                    g_dens(ipm+2), alo, ahi)
            enddo
      enddo !ipm

      stat = .true.
      stat = stat.and.ga_destroy(g_dens_comp_r)
      stat = stat.and.ga_destroy(g_dens_comp_i)
      stat = stat.and.ga_destroy(g_fld1)
      stat = stat.and.ga_destroy(g_fld2)
      stat = stat.and.ga_destroy(g_fldcmplx)
      stat = stat.and.ga_destroy(g_mu_r)
      stat = stat.and.ga_destroy(g_mu_i)
      stat = stat.and.ga_destroy(g_vdimr)
      stat = stat.and.ga_destroy(g_vdimi)
      stat = stat.and.ma_chop_stack(l_dimxyz)
      if (.not.stat) call errquit(
     $      'dimqm_fd: deallocate failed',1,MA_ERR)
      end subroutine dimqm_FD_cycle

      subroutine dim_convert_d2z(n, g_r, g_i, g_z)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"

      integer g_r
      integer g_i
      integer g_z
      integer n, id, i
      double precision valR, valI
      double complex valZ
      id = ga_nodeid()
      if (id.eq.0) then
            do i = 1, n
                  call ga_get(g_r, i, i, 1, 1, valR, 1)
                  call ga_get(g_i, i, i, 1, 1, valI, 1)
                  valZ = dcmplx (valR, valI)
                  call ga_put(g_z, i, i, 1, 1, valZ, 1)
            enddo
      endif
      call ga_sync()
      end subroutine