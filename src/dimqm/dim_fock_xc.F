      subroutine fock_dim(rtdb,geom, nbf_ao, ao_bas_han,
     &   nfock, g_dim, ipm, imag)
c  ==================================================================
c  Wrapper routine to call dim_quadv0 for response calculations
c
c  Called from: dimqm/dimqm_addop.F
c               ddscf/rohf_hessv2.F
c  ==================================================================
      implicit none
#include "errquit.fh"
c
#include "global.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
c
      integer rtdb
      integer geom              ! [input] Geometry handle
      integer nbf_ao,ao_bas_han ! [input] No. of AOs and AO handle
      integer nfock             ! [input] No. of XC matrices
c                               !         = No. of density matx (RDFT)
c                               !         = 4 * No. of density matx (UDFT)
      integer g_dim             ! [input] Array of handles to XC matrices
      integer ipm               ! [input] Control for what dipoles to use
                                !         0 if ignoring +/- directions (static)
                                !         1 for + direction
                                !         2 for - direction
      integer imag              ! [input] More control fro what dipoles to use
                                !         2 when working with imaginary component
                                !         1 when working with real component
c
c     Local declarations
      integer i3
      logical ldebug
      call dimqm_get_debug(rtdb, ldebug)
c
      if (ldebug .and. ga_nodeid().eq.0)
     $   write(luout,*) 'Entered fock_dim'
c
c      dimqm_geom = geom
c      i3 = 3*nDIM
c
c   Call quadrature routine for solving the DIM potential      
c       call dim_grid_quadv0_gen(rtdb, g_dim, nfock,
c     $                          2, .true., ipm, imag)
       call dim_grid_static_gen(rtdb, g_dim)
c
      return
      end subroutine fock_dim
