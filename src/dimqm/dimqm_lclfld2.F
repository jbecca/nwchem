      subroutine dimqm_lclfld2(rtdb,
     $                         g_dipel, 
     $                         g_dipel_i, 
     $                         omega, 
     $                         lifetime)
c       Constructs the local field potential 
c       of the DIM system from the external field
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
#include "cscf.fh"
c-------------------------- Input Variables --------------------------
      integer rtdb
      integer g_dipel
      integer g_dipel_i
      double precision omega
      logical lifetime
c------------------------- Local Variables ----------------------------
      integer id, i3, idir
      integer l_dimxyz, k_dimxyz
      integer l_fld, k_fld
      integer g_fld
      integer g_fld_tot
      logical lfirst
      double complex zsum
      external zsum
      integer dims(3), chunk(3)
      integer g_dim_r, g_dim_i
      integer g_mu_r, g_mu_i
      integer g_temp(2)
      double precision prefactor
      integer alo(3), ahi(3)
      integer blo(2), bhi(2)
      integer clo(3), chi(3)
      logical stat
      
      id = ga_nodeid()
      i3 = 3*DIM_nAtoms
c --------- Allocate arrays and populate ---------
      if(.not.ma_push_get(mt_dbl,i3,'dim coords',l_dimxyz,k_dimxyz))
     $   call errquit('dimqm_lclfld2: malloc k_xyz failed',1,MA_ERR)
      if(.not.rtdb_get(rtdb,'dimqm:coords',mt_dbl,i3,dbl_mb(k_dimxyz)))
     $   call errquit('dimqm_lclfld2: rtdb get failed',1,RTDB_ERR)
      if(.not.ma_push_get(mt_dcpl,i3,'dim fld',l_fld,k_fld))
     $   call errquit('dimqm_lclfld2: malloc k_fld failed',1,MA_ERR)
      if(.not.ga_create(mt_dcpl,i3,1,'g_fld',0,0,g_fld))
     $   call errquit('dimqm_lclfld2: create g_fld failed',0,GA_ERR)
      if(.not.ga_create(mt_dcpl,i3,3,'g_fld_tot',0,0,g_fld_tot))
     $   call errquit('dimqm_lclfld2: create g_fld_tot failed',0,GA_ERR)
      if(.not.ga_create(mt_dbl,i3,3,'g_mu_r',0,0,g_mu_r))
     $   call errquit('dimqm_lclfld2: create g_mu_r failed',0,GA_ERR)
      if(.not.ga_create(mt_dbl,i3,3,'g_mu_i',0,0,g_mu_i))
     $   call errquit('dimqm_lclfld2: create g_mu_i failed',0,GA_ERR)
      call ga_zero(g_fld)
      call ga_zero(g_fld_tot)
      call ga_zero(g_mu_r)
      call ga_zero(g_mu_i)
      lfirst = .true.
      call dimqm_setlfirst(rtdb,lfirst)
      do idir = 1, 3
        call zfill(i3, ZERO_C, dcpl_mb(k_fld),1)
        call zfill(DIM_nAtoms, ONE_C, dcpl_mb(k_fld+idir-1),3)
        call ga_put(g_fld, 1, i3, 1, 1, dcpl_mb(k_fld), 1)
        call dimqm_iterative_FD(rtdb,
     $                          g_fld,
     $                          dbl_mb(k_dimxyz),
     $                          DIM_nAtoms)
        call ga_copy_patch('n',g_fld,     1, i3, 1, 1,
     $                         g_fld_tot, 1, i3, idir, idir)
      enddo

      call convert_z2d(g_fld_tot, g_mu_r, g_mu_i)
      dims(1) = 3
      dims(2) = nbf
      dims(3) = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1

      if (.not.nga_create(mt_dbl, 3, dims, 'lclfld:dim_r+',chunk,
     $      g_dim_r))
     $      call errquit('lclfld: ga_create g_dim_r+ failed',1,GA_ERR)
      if (.not.nga_create(mt_dbl, 3, dims, 'lclfld:dim_i+',chunk,
     $      g_dim_i))
     $      call errquit('lclfld: ga_create g_dim_i+ failed',1,GA_ERR)
      call ga_zero(g_dim_r)
      call ga_zero(g_dim_i)
      call dim_grid_FD_gen(
     $      rtdb,
     $      g_dim_r,
     $      g_mu_r,
     $      1,
     $      1
     $)
      call dim_grid_FD_gen(
     $      rtdb,
     $      g_dim_i,
     $      g_mu_i,
     $      1,
     $      1
     $)
      call ga_symmetrize(g_dim_r)
      call ga_symmetrize(g_dim_i)

      ! recover real and imaginary portion from symmetrization
      dims(1) = nbf
      dims(2) = nbf
      chunk(1) = dims(1)
      chunk(2) = -1

      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      clo(1) = 1
      chi(1) = nbf
      clo(2) = 1
      chi(2) = nbf

      do idir = 1, 3
            alo(1) = idir
            ahi(1) = idir
            clo(3) = idir
            chi(3) = idir
            call nga_add_patch(ONE, g_dipel, clo, chi,
     $                         ONE, g_dim_r, alo, ahi,
     $                              g_dipel, clo, chi)
            call nga_add_patch(ONE, g_dipel_i, clo, chi,
     $                         ONE, g_dim_i,   alo, ahi,
     $                              g_dipel_i, clo, chi)
      enddo

      stat = .true.
      stat = stat.and.ga_destroy(g_mu_r)
      stat = stat.and.ga_destroy(g_mu_i)
      stat = stat.and.ga_destroy(g_dim_r)
      stat = stat.and.ga_destroy(g_dim_i)
      stat = stat.and.ga_destroy(g_fld)
      stat = stat.and.ga_destroy(g_fld_tot)
      stat = stat.and.ma_chop_stack(l_dimxyz)

      end subroutine dimqm_lclfld2