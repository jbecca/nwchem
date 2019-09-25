      subroutine dimqm_addop_cmplx(    rtdb,
     $                                  geom,
     $                                  basis,
     $                                  g_fock, 
     $                                  ncomp, 
     $                                  limag, 
     $                                  lifetime, 
     $                                  nDIM)
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
c
c     Local variables
      integer g_tmp1, g_tmp2, g_dcv
c      integer l_dimxyz, k_dimxyz
      double precision dimxyz(3, nDIM)
c      integer l_muind, k_muind
      double precision muind(3, nDIM, 2)
      integer dims(3), chunk(3)
      character*(255) cstemp
      integer g_pmats(2), g_pmata(2), g_h1mat(2)
      integer ipm
      integer g_dens_r(2)
      integer g_dens_i(2)
      integer alo(3), ahi(3)
      integer blo(2), bhi(2)
      integer g_dens_comp_r
      integer g_dens_comp_i
      integer xend
      double precision pre_factor
      double precision muold(3, nDIM, 2)
      
      double precision dx_r, dy_r, dz_r
      double precision dx_i, dy_i, dz_i
      double precision dsum
      external dsum
      integer i3, ivec, n
c      integer l_fld, k_fld
      double precision fld(3, nDIM, 2)
      integer g_dim_r(2)
      integer g_dim_i(2)
      integer nvir, voff, xoff
      integer nmul, nset, nblock
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      character*(1) direction(3)
      character*(1) pm(2)
      data direction /'x', 'y', 'z'/
      data pm /'+', '-'/
      logical firsttime
c      double precision screen(nDIM)
      double precision dimErr(3,2,2)
      double precision calcErr
      external calcErr
      integer id
c
      id = ga_nodeid()
      i3 = nDIM * 3
C      g_tmp1 = ga_create_atom_blocked(geom, basis, 'dim_addop: tmp1')
c
c   Allocate new arrays
c      if(.not.ma_push_get(mt_dbl,i3*2,'dimrsp:fld',l_fld,k_fld))
c     $  call errquit('malloc dimrsp:fld failed',1,MA_ERR)
c
c      if(.not.ma_push_get(mt_dbl,i3*2,'dimrsp:muind',
c     $                                            l_muind,k_muind))
c     $  call errquit('malloc dimrsp:muind failed',1,MA_ERR)
c
c      if(.not.ma_push_get(mt_dbl,i3,'dimrsp:xyz',l_dimxyz,k_dimxyz))
c     $  call errquit('malloc dimrsp:xyz failed',1,MA_ERR)
c
      if(.not.rtdb_get(rtdb,'dimpar:coords', mt_dbl, i3, dimxyz))
     $  call errquit('get dimpar:coords failed', 1, RTDB_ERR)
c
      g_dens_comp_r = ga_create_atom_blocked(geom,basis,
     $                          'real density matrix comp')
      if (lifetime) then
        g_dens_comp_i = ga_create_atom_blocked(geom,basis,
     $                          'imag density matrix comp')
      end if
cc
cc   =============================
cc   Solve for induced dipoles +/-
cc   =============================
cc
cc     Loop over perturbations
c      do n = 1, 3
c        do ipm = 1, 2
c          call ga_zero(g_dens_comp_r)
c          if (lifetime) call ga_zero(g_dens_comp_i)
c          alo(1) = n
c          ahi(1) = n
cc
cc       Copy current perturbation into g_dens_comp
c          call nga_copy_patch('N',g_dens_r(ipm), alo, ahi,
c     $                          g_dens_comp_r, blo, bhi)
c          if (lifetime) then
c            call nga_copy_patch('N',g_dens_i(ipm), alo, ahi,
c     $                              g_dens_comp_i, blo, bhi)
c          end if
c          muind = ZERO
c          fld = ZERO
c          firsttime = .false.
c          if(.not.rtdb_get(rtdb,
c     $                'dimqm:muind_'//direction(n)//'_r'//pm(ipm),
c     $                              mt_dbl, i3, muold(:,:,1))) then
c            if(id.eq.0) write(luout,*) "First cycle, no old dipoles!"
c            muold = ZERO
c            firsttime = .true.
c            dimqm_seeded = .false.
cc            xyz_seeded(3*(n-1)+ipm) = .false.
c            if(dimtol0 < 1.0d-4 .and. .not. dimqm_noseed) then
c              dimtolxyz(ipm*3 - 1 + n) = 1.0d-4
c              if(id.eq.0) then
c                write(luout,*) "Requested tolerance below 1.0d-4"
c                write(luout,*) "Setting "//direction(n)//pm(ipm)//
c     $                         " dir tolerance to 1.0d-4 to start"
c              end if
c            end if
c          else
c            if(.not.rtdb_get(rtdb,
c     $                'dimqm:muind_'//direction(n)//'_i'//pm(ipm),
c     $                           mt_dbl, i3, muold(:,:,2)))
c     $          call errquit('get dimqm:muold failed',1,RTDB_ERR)
c          end if
cc         Set convergence tolerance
cc          dimtol = dimtolxyz(ipm*3 - 1 + n)
cc          dimqm_seeded = xyz_seeded(ipm*3 - 1 + n)
cc          dimtol = 1.0d-7
cc          dimqm_noseed = .true.
cc          call dfill(i3*2, ZERO, dbl_mb(k_muind), 1)
cc          call dfill(i3*2, ZERO, dbl_mb(k_fld), 1)
cc
cc       Real portion of E-Field
cc        write(luout,*) "REAL"
cc        call ga_print(g_dens_comp_r)
c          call dimqm_EqmE(rtdb, g_dens_comp_r, geom, basis,
c     $               fld(:,:,1), dimxyz)
cc
cc       Imaginary portion of E-Field
cc        write(luout,*) "IMAG"
cc        call ga_print(g_dens_comp_i)
c          if (lifetime) then
c            call dimqm_EqmE(rtdb, g_dens_comp_i, geom, basis,
c     $                      fld(:,:,2), dimxyz)
c          end if
cc
cc       Solve for induced dipoles
c          call dimqm_f2d(rtdb, fld, muind, muold, dimxyz, 2,
c     $                   direction(n), pm(ipm))
cc
cc         Write induced dipoles to RTDB
c          dx_r = SUM(muind(1,:,1))
c          dy_r = SUM(muind(2,:,1))
c          dz_r = SUM(muind(3,:,1))
c          dx_i = SUM(muind(1,:,2))
c          dy_i = SUM(muind(2,:,2))
c          dz_i = SUM(muind(3,:,2))
c          if(id.eq.0) then
c            write(luout,*) "Total induced dipole moment for "//
c     $                  direction(n)//pm(ipm)//" perturbation"
c            write(luout,*) "X:", dx_r, dx_i
c            write(luout,*) "Y:", dy_r, dy_i
c            write(luout,*) "Z:", dz_r, dz_i
c            write(luout,*) ''
c          end if
c          dimErr(n, ipm, 1) = calcErr(i3, muold(:,:,1), muind(:,:,1))
c          dimErr(n, ipm, 2) = calcErr(i3, muold(:,:,2), muind(:,:,2))
c          if(id.eq.0) then
c            write(luout,*) "Max error in real dipoles:",
c     $                       dimErr(n, ipm, 1)
c            write(luout,*) "Max error in imag dipoles:",
c     $                       dimErr(n, ipm, 2)
c          end if
cc          if(dimErr(n, ipm, 1)/dimtol < HUNDRED
cc     $              .and. dimErr(n, ipm, 2)/dimtol < HUNDRED
cc     $              .and. .not. xyz_seeded(ipm*3 - 1 + n)
cc     $              .and. .not. firsttime) then
cc            xyz_seeded(ipm*3 - 1 + n) = .true.
cc            write(luout,*) "Error within 10^2 of", dimtol, "for "//
cc     $                     direction(n)//pm(ipm)//" dir"
cc            write(luout,*) "Setting current "//direction(n)//pm(ipm)//
cc     $                     " dir as seed"
cc            write(luout,*)"Reverting tolerance back to", dimtol0
cc            dimtolxyz(ipm*3 - 1 + n) = dimtol0
cc          end if
c          if(.not.rtdb_put(rtdb,
c     $                'dimqm:muind_'//direction(n)//'_r'//pm(ipm),
c     $                              mt_dbl, i3, muind(:,:,1)))
c     $        call errquit('put dimqm:muind_p failed',1,RTDB_ERR)
c          if(.not.rtdb_put(rtdb,
c     $                'dimqm:muind_'//direction(n)//'_i'//pm(ipm),
c     $                           mt_dbl, i3, muind(:,:,2)))
c     $        call errquit('put dimqm:muind_p failed',1,RTDB_ERR)
c        end do ! ipm = 1, 2
c      end do ! ivec = 1, 3
cc      if(MAXVAL(dimErr) <= 1.0d-4) then
cc        write(luout,*) "Dipole error below 1d-4"
cc        write(luout,*) "Shutting down DIM"
cc        dimqm_on = .false.
cc      end if
cc
cc   Destroy GAs we don't need anymore
c      if (.not. ga_destroy(g_dens_comp_r)) call errquit
c     $    ('addop: dens_comp_r GA?',0, GA_ERR)
c      if (lifetime) then
c      if (.not. ga_destroy(g_dens_comp_i)) call errquit
c     $    ('addop: dens_comp_i GA?',0, GA_ERR)
c      end if
c      do ipm = 1,2
c        if (.not. ga_destroy(g_dens_r(ipm))) call errquit
c     $     ('addop: dens_r GA?',0, GA_ERR)
c        if (lifetime) then
c        if (.not. ga_destroy(g_dens_i(ipm))) call errquit
c     $     ('addop: dens_i GA?',0, GA_ERR)
c        end if
c      end do
cc
cc   Deallocate l_fld, l_muind, l_dimxyz
cc      if (.not. ma_chop_stack(l_fld)) call errquit
cc     $   ('addop: fld MA?', 0, MA_ERR)
cc
cc   ====================================================
cc   Solve for DIM potential, both real and imaginary S/A
cc   ====================================================
cc
c      dims(1) = 3
c      dims(2) = nbf
c      dims(3) = nbf
c      chunk(1) = dims(1)
c      chunk(2) = -1
c      chunk(3) = -1
cc
cc   Real +
c      if (.not. nga_create (MT_DBL, 3, dims, 'addop:dim_r+',chunk,
c     &                        g_dim_r(1)))
c     &   call errquit('addop: could not allocate g_dim_r+',1,GA_ERR)
c      call ga_zero(g_dim_r(1)) 
c      call fock_dim(geom, nbf, basis, 3, g_dim_r(1), 1, 1)
c      call ga_symmetrize(g_dim_r(1))
cc
cc   Real -
c      if (.not. nga_create (MT_DBL, 3, dims, 'addop:dim_r-',chunk,
c     &                        g_dim_r(2)))
c     &   call errquit('addop: could not allocate g_dim_r-',1,GA_ERR)
c      call ga_zero(g_dim_r(2))
c      call fock_dim(geom, nbf, basis, 3, g_dim_r(2), 2, 1)
c      call ga_antisymmetrize(g_dim_r(2))
c      if (lifetime) then
cc
cc   Imaginary +
c      if (.not. nga_create (MT_DBL, 3, dims, 'addop:dim_i+',chunk,
c     &                        g_dim_i(1)))
c     &   call errquit('addop: could not allocate g_dim_i+',1,GA_ERR)
c      call ga_zero(g_dim_i(1))
c      call fock_dim(geom, nbf, basis, 3, g_dim_i(1), 1, 2)
c      call ga_symmetrize(g_dim_i(1))
cc
cc   Imaginary -
c      if (.not. nga_create (MT_DBL, 3, dims, 'addop:dim_i-',chunk,
c     &                        g_dim_i(2)))
c     &   call errquit('addop: could not allocate g_dim_i-',1,GA_ERR)
c      call ga_zero(g_dim_i(2))
c      call fock_dim(geom, nbf, basis, 3, g_dim_i(2), 2, 2)
c      call ga_antisymmetrize(g_dim_i(2))
c      end if
cc
cc   ======================================
cc   Undo the symmetrization to recover +/-
cc   ======================================
cc
c      do ivec = 1, 3
c        alo(1) = ivec
c        ahi(1) = ivec
cc       ************
cc       Real portion
cc       ************
c        call nga_copy_patch ('N',g_dim_r(1),alo,ahi,g_pmats(1),blo,bhi)
c        call nga_copy_patch ('N',g_dim_r(2),alo,ahi,g_pmats(2),blo,bhi)
cc
cc       it might be necessary to use 0.5 here instead of 1.0
cc       (note: that turned out NOT to be the case after some testing)
c        pre_factor = 1.0d0
c        call ga_sync()
c        if (.not.limag) then
cc         real perturbation: 
c          call nga_add_patch (pre_factor, g_pmats(1), blo, bhi,
c     &       pre_factor, g_pmats(2), blo, bhi,
c     &       g_dim_r(1), alo, ahi)
c          call nga_add_patch (pre_factor, g_pmats(1), blo, bhi,
c     &       -pre_factor, g_pmats(2), blo, bhi,
c     &       g_dim_r(2), alo, ahi)
c        else
cc         imaginary perturbation:
c          call nga_add_patch (pre_factor, g_pmats(1), blo, bhi,
c     &       pre_factor, g_pmats(2), blo, bhi,
c     &       g_dim_r(1), alo, ahi)
c          call nga_add_patch (-pre_factor, g_pmats(1), blo, bhi,
c     &       pre_factor, g_pmats(2), blo, bhi,
c     &       g_dim_r(2), alo, ahi)
c        end if  ! if .not.limag
c        if (lifetime) then
cc       *****************
cc       Imaginary portion
cc       *****************
c        call nga_copy_patch ('N',g_dim_i(1),alo,ahi,g_pmats(1),blo,bhi)
c        call nga_copy_patch ('N',g_dim_i(2),alo,ahi,g_pmats(2),blo,bhi)
cc
cc       it might be necessary to use 0.5 here instead of 1.0
cc       (note: that turned out NOT to be the case after some testing)
c        pre_factor = 1.0d0
c        call ga_sync()
c        if (.not.limag) then
cc         real perturbation: 
c          call nga_add_patch (pre_factor, g_pmats(1), blo, bhi,
c     &       pre_factor, g_pmats(2), blo, bhi,
c     &       g_dim_i(1), alo, ahi)
c          call nga_add_patch (pre_factor, g_pmats(1), blo, bhi,
c     &       -pre_factor, g_pmats(2), blo, bhi,
c     &       g_dim_i(2), alo, ahi)
c        else
cc         imaginary perturbation:
c          call nga_add_patch (pre_factor, g_pmats(1), blo, bhi,
c     &       pre_factor, g_pmats(2), blo, bhi,
c     &       g_dim_i(1), alo, ahi)
c          call nga_add_patch (-pre_factor, g_pmats(1), blo, bhi,
c     &       pre_factor, g_pmats(2), blo, bhi,
c     &       g_dim_i(2), alo, ahi)
c        end if  ! if .not.limag
c        end if ! lifetime
c      enddo                     ! ivec = 1,nvec
c100   continue
cc
cc   ====================================
cc   Add DIM potential to the Fock matrix
cc   ====================================
cc
c      g_dcv = ga_create_atom_blocked(geom, basis, 'rohf_h2e3: dcv')
c      xoff = 1
c      voff = nclosed + nopen + 1
c      xend = nvir*nclosed
c      do ivec = 1, 3 ! Loop over perturbations
c        alo(1) = ivec
c        ahi(1) = ivec
c        do ipm = 1, 2! Loop over +/-
cc         We only add the + direction of the DIM potential to both +/- of the Fock matrix
cc   Real Portion
c          call nga_copy_patch('N',g_dim_r(ipm),alo,ahi,g_dcv,blo,bhi)
c          call ga_scale(g_dcv, two)
c          call ga_matmul_patch('n', 'n', two, zero,
c     $                           g_dcv,   1, nbf, 1, nbf,
c     $                           g_movecs, 1, nbf, 1, nclosed,
c     $                           g_tmp1,  1, nbf, 1, nclosed)
c          call ga_sync()
c          call ga_matmul_patch('t', 'n', one, zero,
c     $                           g_movecs, voff, nmo, 1, nbf,
c     $                           g_tmp1, 1, nbf,  1, nclosed,
c     $                           g_tmp2, 1, nvir, 1, nclosed)
c          call ga_sync()
c          call ga_mat_to_vec(g_tmp2, 1, nvir, 1, nclosed, g_Ax_r(ipm),
c     $                         xoff, ivec, four, '+')
cc
cc   Imaginary Portion
c          if (lifetime) then
c          call nga_copy_patch('N',g_dim_i(ipm),alo,ahi,g_dcv,blo,bhi)
c          call ga_scale(g_dcv, two)
c          call ga_matmul_patch('n', 'n', two, zero,
c     $                           g_dcv,   1, nbf, 1, nbf,
c     $                           g_movecs, 1, nbf, 1, nclosed,
c     $                           g_tmp1,  1, nbf, 1, nclosed)
c          call ga_sync()
c          call ga_matmul_patch('t', 'n', one, zero,
c     $                           g_movecs, voff, nmo, 1, nbf,
c     $                           g_tmp1, 1, nbf,  1, nclosed,
c     $                           g_tmp2, 1, nvir, 1, nclosed)
c          call ga_sync()
cc
cc   ******NOTE JEM******
cc   We can remove the - sign if we apply the sign change discussed in aoresponse and rohf_hessv_xx3
cc   ********************
cc
c          call ga_mat_to_vec(g_tmp2, 1, nvir, 1, nclosed, g_Ax_i(ipm),
c     $                         xoff, ivec, four, '+')
c          end if !lifetime
c        end do !ipm = 1, 2
c      end do !ivec = 1, 3
cc   ========
cc   Clean up
cc   =======
c      lfirst = .false.
c      do ipm = 1,2
c        if (.not. ga_destroy(g_pmats(ipm))) call errquit
c     $     ('addop: pmats GA?', 0, GA_ERR)
c        if (.not. ga_destroy(g_pmata(ipm))) call errquit
c     $     ('addop: pmata GA?', 0, GA_ERR)
c        if (.not. ga_destroy(g_h1mat(ipm))) call errquit
c     $     ('addop: h1mat GA?', 0, GA_ERR)
c        if (.not. ga_destroy(g_dim_r(ipm))) call errquit
c     $     ('addop: dim_r GA?', 0, GA_ERR)
c        if (lifetime) then
c        if (.not. ga_destroy(g_dim_i(ipm))) call errquit
c     $     ('addop: dim_i GA?', 0, GA_ERR)
c        end if
c      enddo                     ! ipm = 1,2
cc
c      if (.not. ga_destroy(g_tmp1)) call errquit
c     $   ('addop: tmp1 GA?', 0, GA_ERR)
c      if (.not. ga_destroy(g_tmp2)) call errquit
c     $   ('addop: tmp2 GA?', 0, GA_ERR)
c      if (.not. ga_destroy(g_dcv)) call errquit
c     $   ('addop: dcv GA?',0, GA_ERR)
cc
      end subroutine dimqm_addop_cmplx
