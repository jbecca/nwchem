      subroutine dimqm_FD_cycle(    rtdb,
     $                              geom,
     $                              basis,
     $                              g_fock, 
     $                              ncomp, 
     $                              limag, 
     $                              lifetime, 
     $                              nDIM,
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
      integer rtdb, geom, basis
      integer g_fock(ncomp)   ! A matrix handle (complex) [IN]
      integer ncomp
      logical limag        ! Imaginary perturbation?   [IN]
      logical lifetime     ! Damping or no damping     [IN]
      integer nDIM         ! number of DIM atoms [IN]
      integer g_dens(4)
c
      integer i3
      integer id, n, ipm
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      integer g_dens_comp_i, g_dens_comp_r
c
      id = ga_nodeid()
      i3 = DIM_nAtoms * 3

      g_dens_comp_r = ga_create_atom_blocked(geom,basis,
     $                          'real density matrix comp')
      if (lifetime) then
        g_dens_comp_i = ga_create_atom_blocked(geom,basis,
     $                          'imag density matrix comp')
      end if

      do n = 1, 3 ! loop over directions
            do ipm = 1, 2 ! loop over perturbation sign
            enddo !ipm 
      enddo !n
      !     call EqmE for real part of field
      !      call dimqm_eqme
      !     call EqmE for imag part of field
      !      call dimqm_eqme
      !     call f2d for real and imag induced dipoles
      !      call dimqm_f2d

      if (.not. ga_destroy(g_dens_comp_r)) call errquit
     $    ('addop: dens_comp_r GA?',0, GA_ERR)
      if (lifetime) then
      if (.not. ga_destroy(g_dens_comp_i)) call errquit
     $    ('addop: dens_comp_i GA?',0, GA_ERR)
      end if

      end subroutine dimqm_FD_cycle
