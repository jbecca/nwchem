      subroutine dimqm_pol_static(rtdb)
c  ==================================================================
c  Calculate the standalone polarizability of the DIM system.  Here,
c  we apply a unit field in each cartesian direction and determine
c  the DIM system's response.
c
c  Called from: property/aoresponse_driver.F
c  ==================================================================
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
#include "util.fh"
#include "dimqm_constants.fh"
#include "dimqm.fh"
c
c     Input variables
      integer rtdb        ! RTDB handle
c
c     MA handles
      integer k_fld, l_fld
      integer l_dimxyz, k_dimxyz ! DIM coordinates handles
      integer l_muind,  k_muind  ! DIM dipoles handles
      integer l_muold,  k_muold
c
c     Local variables
      integer idir, m, i3
      integer icmplx
      double precision pol(3,3,2)
      double precision time
      logical stat
      integer id
      double precision dsum
      external dsum
      logical seed_save
      logical dimqm_noseed
      integer nDIM, nOrder
      logical lfirst
      logical ldebug
      integer g_fld
c
c   Common variables used
c      logical lcpim, lcmplx, ldebug
c      integer nOrder, nDIM
c
c   ==========
c   Initialize 
c   ==========
c
      call dimqm_getnDIM(rtdb, nDIM)
      i3 = 3*DIM_nAtoms
      call dimqm_get_debug(rtdb, ldebug)
      id = ga_nodeid()
      if(id.eq.0.and.ldebug) write(LuOut,*)
     $   "Begin dimqm_pol_static routine"
c
c   ========================
c   Allocate arrays and fill
c   ========================
c
c     Coords
      if(.not.ma_push_get(mt_dbl,i3,'dim coords',l_dimxyz,k_dimxyz))
     $   call errquit('dimqm_pol malloc k_xyz failed',1,MA_ERR)
      if(.not.rtdb_get(rtdb,'dimqm:coords',mt_dbl,i3,dbl_mb(k_dimxyz)))
     $   call errquit('dimqm:coords rtdb get failed',1,RTDB_ERR)
c
c     GAs
      if(.not.ga_create(mt_dbl,DIM_nOrder,1,'g_fld',0,0,g_fld))
     $   call errquit('dimqm_gs: create g_fld failed',0,GA_ERR)
      call ga_zero(g_fld)
c     Induced Dipoles
c      if(.not.ma_push_get(mt_dbl,i3*icmplx,'dim muind',
c     $                                            l_muind,k_muind))
c     $   call errquit('dimqm_pol malloc k_muind failed',1,MA_ERR)
c      if(.not.ma_push_get(mt_dbl,i3*icmplx,'dim:muold',
c     $                                            l_muold,k_muold))
c     $   call errquit('dimqm_pol malloc k_muind failed',1,MA_ERR)
c
c     Applied Field
      if(.not.ma_push_get(mt_dbl,i3,'dim fld',l_fld,k_fld))
     $   call errquit('dimqm_pol malloc k_fld failed',1,MA_ERR)
c
c   Zero polarizability
      pol = ZERO
c
c   ===================================================
c   Loop over each unit field and sum up the dipoles to
c   get the polarizability
c   ===================================================
c
      time = util_timer()
      lfirst = .true.
      call dimqm_setlfirst(rtdb, lfirst)
      seed_save = dimqm_noseed
      dimqm_noseed = .true.
      if(.not.rtdb_put(rtdb,'dimqm:noseed',mt_log,1,dimqm_noseed))
     $  call errquit('dimqm_input: noseed put failed', 1, RTDB_ERR)
      do idir = 1, 3
c
c           Zero field
          call dfill(i3, ZERO, dbl_mb(k_fld), 1)
c
c           Unit field in specified direction
          call dfill(DIM_nAtoms, 1.0d0, dbl_mb(k_fld+idir-1), 3)
          call ga_put(g_fld, 1, DIM_nOrder, 1, 1, dbl_mb(k_fld), 1)
          call ga_sync()

          call dimqm_iterative2(rtdb, g_fld, dbl_mb(k_dimxyz),
     $             DIM_nAtoms)
          call ga_get(g_fld, 1, DIM_nOrder, 1, 1, dbl_mb(k_fld), 1)
c           Sum up dipoles to get polarizability
          pol(1,idir,1) = dsum(DIM_nAtoms, dbl_mb(k_fld), 3)
          pol(2,idir,1) = dsum(DIM_nAtoms, dbl_mb(k_fld+1), 3)
          pol(3,idir,1) = dsum(DIM_nAtoms, dbl_mb(k_fld+2), 3)
      end do
c
c   Print polarizability
      if(id.eq.0) then
          write(luout,5000) "Static",
     $                      pol(1,1,1), pol(1,2,1), pol(1,3,1),
     $                      pol(2,1,1), pol(2,2,1), pol(2,3,1),
     $                      pol(3,1,1), pol(3,2,1), pol(3,3,1)
          write(luout,*)
        call util_flush(LuOut)
      end if
      dimqm_noseed = seed_save
c
c   ==============================
c   Deallocate arrays and clean up
c   ==============================
c
      stat = ma_chop_stack(l_dimxyz)
c  jbecca START: dont think we need to create and destroy int mat here.
c        this is called from ao_driver_new, where it is already made
c      call dimqm_tidy()
c  jbecca END
      if(id .eq. 0) write(LuOut,*)
     $  "Total time:", util_timer() - time
 5000 format("Polarizability tensor for DIM system: ", A ,/
     $       1x,'             X              Y              Z',/
     $       1x,'-----------------------------------------------',/
     $       1x,'X ',3f15.7,/
     $       1x,'Y ',3f15.7,/
     $       1x,'Z ',3f15.7,/
     $       1x,'-----------------------------------------------')


      end subroutine dimqm_pol_static

