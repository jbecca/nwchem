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
c------------------------- Local Variables ----------------------------
      integer id, i3, idir
      integer l_dimxyz, k_dimxyz
      integer l_fld, k_fld
      integer g_fld
      integer g_fld_tot
      logical lfirst
      
      id = ga_nodeid()
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
      call ga_zero(g_fld)
      call ga_zero(g_fld_tot)
      lfirst = .true.
      call dimqm_setlfirst(rtdb,lfirst)
      do idir = 1, 3
        call zfill(i3, ZERO_C, dcpl_mb(k_fld),1)
        call zfill(DIM_nAtoms, ONE_C, dcpl_mb(k_fld+idir-1),3)
        call ga_put(g_fld, 1, i3, 1, 1, dcpl_mb(kfld), 1)
        call dimqm_iterative_FD(rtdb,
     $                          g_fld,
     $                          dbl_mb(k_fld),
     $                          DIM_nAtoms)
        call nga_copy_patch('n' g_fld,     1, i3, 1, 1,
     $                          g_fld_tot, 1, i3, idir, idir)
        call ga_print(g_fld)
      enddo
      call ga_print(g_fld_tot)
      end subroutine dimqm_lclfld2