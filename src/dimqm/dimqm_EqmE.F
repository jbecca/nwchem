      subroutine dimqm_EqmE(rtdb, g_dens, geom, basis,
     $                      fld, xyz, nDIM)
c  
c      Author: Justin Moore
c 
c      Called from: dimqm_main.F dimqm_addop.F dimqm_indDipoles.F
c                   dimqm_Dens2Dipoles.F
c
c      Calls: hnd_elfcon.F
c
c      Subroutine to calculate the electric field due to the QM
c      electrons
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
#include "geom.fh"
c    Input Variables
      integer rtdb
      integer g_dens
      integer geom
      integer basis
      double precision fld(3, nDIM)
      double precision xyz(3, nDIM)
      integer nDIM
c
      double precision scrn(nDIM)
      integer i, j, n, id
      logical stat
      integer nNuc
      integer l_qmxyz, k_qmxyz
      integer l_dimxyz, k_dimxyz
      integer l_dimfld, k_dimfld
c      double precision dimfactors(nDIM)
      double precision temp
      double precision pref
      double precision ox, oy, oz
      double precision ac, aq, aa
      logical ldebug, status

c      ac = 1.0d0
c      aq = 1.0d0
c      aa = (ac*aq)/(ac + aq)
c      pref = sqrt(PI/(ac+aq))
      status = rtdb_parallel(.true.)
      id = ga_nodeid()
      call dimqm_get_debug(rtdb,ldebug)
      if(ldebug) then
        write(luout,*) "Entering EqmE"
      end if
c      if(.not. geom_ncent(geom, nNuc))
c     $  call errquit("EqmE: geom_ncent", 1, GEOM_ERR)
c      if(ldebug) write(Luout,*)'got QM ncent'
c      if(.not. ma_push_get(mt_dbl,3*nNuc,'QM xyz',l_qmxyz,k_qmxyz))
c     $  call errquit("EqmE: malloc k_qmxyz", 1, MA_ERR)
      if(.not. ma_push_get(mt_dbl,3*nDIM,'DIM xyz',l_dimxyz,k_dimxyz))
     $  call errquit("EqmE: malloc k_dimxyz", 1, MA_ERR)
      do i = 1, nDIM
        dbl_mb(k_dimxyz+3*(i-1)+0) = xyz(1,i)
        dbl_mb(k_dimxyz+3*(i-1)+1) = xyz(2,i)
        dbl_mb(k_dimxyz+3*(i-1)+2) = xyz(3,i)
      enddo
c  try allocating field to see if this fixes the issue for hnd_elfcon
      if(.not. ma_push_get(mt_dbl,3*nDIM,'DIM fld',l_dimfld,k_dimfld))
     $  call errquit("EqmE: malloc k_dimfld", 1, MA_ERR)

      if(ldebug) write(Luout,*)'copied coords to k_dimxyz'
c      if(ldebug) write(Luout,*)'allocated QM geom'
c      if(.not. geom_cart_coords_get(geom, dbl_mb(k_qmxyz)))
c     $  call errquit("EqmE: qm coords get failed", 1, GEOM_ERR)
c      if(ldebug) write(Luout,*)'got QM coords'
      call ga_sync()
      call hnd_elfcon(basis, geom, g_dens, dbl_mb(k_dimxyz),
     $      nDIM,dbl_mb(k_dimfld),1)
      call ga_sync()
c      call hnd_elfcon(basis, geom, g_dens, xyz,
c     $      nDIM,fld,1)
      if(ldebug) write(Luout,*)'back from hnd_elfcon'
c      call ga_sync()
      if(ldebug) write(Luout,*)'ga have been synced'
c      if(.not.ma_pop_stack(l_qmxyz))
c     %  call errquit('problem destroying qmxyz in EqmE', 1, MA_ERR)
c      if(.not.ma_pop_stack(l_dimxyz))
c     %  call errquit('problem destroying dimxyz in EqmE', 1, MA_ERR)
      if(.not.ma_pop_stack(l_dimfld))
     %  call errquit('problem destroying dimfld in EqmE', 1, MA_ERR)
      if(ldebug) then
        write(luout,*) "Exiting EqmE"
      end if
      end subroutine dimqm_EqmE
