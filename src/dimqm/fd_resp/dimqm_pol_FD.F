      subroutine dimqm_pol_FD(rtdb, om)
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
      double precision om ! Frequency of incident light (zero if static calculation)
c
c     MA handles
      integer k_fld, l_fld
      integer l_dimxyz, k_dimxyz ! DIM coordinates handles
c
c     Local variables
      integer idir, i3, j
      double precision pol(3,3,2)
      double complex polc(3,3)
      double precision time
      logical stat
      integer id
      double precision dsum
      external dsum
      double complex zsum
      external zsum
      logical lfirst
      logical ldebug
c     GAs
      integer g_fld
c     Init
      i3 = 3 * DIM_nAtoms
      call dimqm_get_debug(rtdb, ldebug)
      id = ga_nodeid()

c     Allocate arrays and fill
      if(.not.ma_push_get(mt_dbl,i3,'dim coords',l_dimxyz,k_dimxyz))
     $  call errquit('dimqm_fd_pol malloc k_xyz failed',1,MA_ERR)
      if(.not.rtdb_get(rtdb,'dimqm:coords',mt_dbl,i3,dbl_mb(k_dimxyz)))
     $  call errquit('dimqm:coords rtdb get failed',1,RTDB_ERR)

      if(.not.ma_push_get(mt_dcpl,i3,'dim fld',l_fld,k_fld))
     $  call errquit('dimqm_pol_FD malloc k_fld failed',1,MA_ERR)

      if(.not.ga_create(mt_dcpl,i3,1,'g_fld',0,0,g_fld))
     $      call errquit('dimqm_fd_pol: create g_fld failed',0,GA_ERR)

      call ga_zero(g_fld)

      pol = ZERO

      time = util_timer()
      lfirst = .true.
      call dimqm_setlfirst(rtdb, lfirst)

      do idir = 1, 3
            call zfill(i3, ZERO_C, dcpl_mb(k_fld), 1)
            call zfill(DIM_nAtoms, ONE_C, dcpl_mb(k_fld+idir-1), 3)
            call ga_put(g_fld, 1, i3, 1, 1, dcpl_mb(k_fld), 1)
            call ga_sync()
            call dimqm_iterative_FD(rtdb,
     $                              g_fld,
     $                              dbl_mb(k_dimxyz),
     $                              DIM_nAtoms)
            call ga_get(g_fld, 1, i3, 1, 1, dcpl_mb(k_fld), 1)
            polc(1,idir) = zsum(DIM_nAtoms, dcpl_mb(k_fld),3)
            polc(2,idir) = zsum(DIM_nAtoms, dcpl_mb(k_fld+1),3)
            polc(3,idir) = zsum(DIM_nAtoms, dcpl_mb(k_fld+2),3)
      enddo

      if(id.eq.0) then
            write(luout,5000) "Frequency-Dependent Real",
     $         REAL(polc(1,1)), REAL(polc(1,2)), REAL(polc(1,3)),
     $         REAL(polc(2,1)), REAL(polc(2,2)), REAL(polc(2,3)),
     $         REAL(polc(3,1)), REAL(polc(3,2)), REAL(polc(3,3))
            write(luout,*)
            write(luout,5000) "Frequency-Dependent Imag",
     $         AIMAG(polc(1,1)), AIMAG(polc(1,2)), AIMAG(polc(1,3)),
     $         AIMAG(polc(2,1)), AIMAG(polc(2,2)), AIMAG(polc(2,3)),
     $         AIMAG(polc(3,1)), AIMAG(polc(3,2)), AIMAG(polc(3,3))
            write(luout,*)
            call util_flush(LuOut)
      endif

      if (.not.ga_destroy(g_fld))
     $   call errquit('ga destroy g_fld failed',0,GA_ERR)
 5000 format("Polarizability tensor for the DIM system: ", A ,/
     $      1x,'              X                 Y                 Z',/
     $      1x,'----------------------------------------------------',/
     $      1x,'X ',3f15.7,/
     $      1x,'Y ',3f15.7,/
     $      1x,'Z ',3f15.7,/
     $      1x,'----------------------------------------------------')

      end subroutine dimqm_pol_FD