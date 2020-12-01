      subroutine dimqm_addop2(rtdb, geom, basis, nbf,
     $                        nmo, nclosed, nopen, g_movecs,
     $                        g_x_r, g_x_i, g_Ax_r, g_Ax_i,
     $                       ncomp, limag, lifetime, nDIM)
c
c     Author: Jeff Becca
c
c     Called from: rohf_hessv3.F rohf_hessv3_ext.F 
c
c     Subroutines called from: CalcPerturbedTDPmat1_opt.F, dimqm_EqmE.F
c                              dimqm_f2d.F dim_fock_xc.F
c   
c     Calculates and adds the frequency-dependent DIM potential to the 
c     response Fock matricies (real and imaginary).  Requires knowledge
c     of both the real and imaginary vectors simultaneously which is
c     why this has to be located here, unlike the static routine.
c
c    It seems that a speed up opportunity is present by passing in the 
c    perturbed density that is calculated in some previous routines.
c    If this can be done, can avoid a large part of the start of this 
c    routine that gets the symmetric and antisymmetric real and imaginary
c    perturbed densities
c
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "dimqm_constants.fh"
#include "dimqm.fh"
c
c     Input Variables
      integer rtdb
      integer g_x_r(2)     ! A matrix handle (real)      [IN]
      integer g_x_i(2)     ! A matrix handle (imaginary) [IN]
      integer g_Ax_r(2)    ! F matrix handle (real)      [IN/OUT]
      integer g_Ax_i(2)    ! F matrix handle (imaginary) [IN/OUT]
      integer ncomp        ! num of components (+/-)     [IN]
      logical limag        ! Imaginary perturbation?     [IN]
      logical lifetime     ! Damping or no damping
      integer nDIM
      integer geom, basis, nbf
      integer nmo, nclosed, nopen, g_movecs
c
c     Local variables
      integer g_tmp1, g_tmp2, g_dcv
      logical stat
      double precision dimxyz(3, nDIM)
      integer l_dimxyz, k_dimxyz
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
      integer g_fld1, g_fld2, g_fldcmplx
      integer g_fld_tot, g_mu_r(2), g_mu_i(2)
      double precision pre_factor
      double precision muold(3, nDIM, 2)
      
      double precision dx_r, dy_r, dz_r
      double precision dx_i, dy_i, dz_i
      double precision dsum, rmax
      external dsum
      integer i3, ivec, n
c      integer l_fld, k_fld
      double precision fld(3, nDIM, 2)
      integer g_dim_r(2)
      integer g_dim_i(2)
      integer nvir, voff, xoff
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      character*(1) direction(3)
      character*(1) pm(2)
      data direction /'x', 'y', 'z'/
      data pm /'+', '-'/
      logical firsttime
c      double precision screen(nDIM)
      double precision dimErr(3,2,2)
      double precision calcErr, dimtol, dimtol0, dimtolxyz(6)
      external calcErr
      integer id
      logical dimqm_noseed, dimqm_seeded
      logical lfirst
      logical ldebug
c
      id = ga_nodeid()
      call dimqm_get_debug(rtdb, ldebug)
      if (ldebug .and. id .eq. 0) then
        write(luout,*) "Start dimqm_addop"
      end if
      call dimqm_getlfirst(rtdb, lfirst)
      nvir = nmo - nclosed - nopen
      i3 = nDIM * 3
      g_tmp1 = ga_create_atom_blocked(geom, basis, 'dim_addop: tmp1')
      g_tmp2 = ga_create_atom_blocked(geom, basis, 'dim_addop: tmp2')

      dims(1) = nbf
      dims(2) = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
c
c   Allocate GAs
      if (ldebug .and. id .eq. 0) then
        write(luout,*) "Start allocation"
      end if
      do ipm = 1, 2
        write(cstemp,'(a,i1)') 'pmats_',ipm
        if (.not.nga_create(MT_DBL,2,dims,cstemp(1:7),chunk,
     &     g_pmats(ipm))) call
     &     errquit('dim_addop: nga_create failed '//cstemp(1:7),
     &     0,GA_ERR)
        call ga_zero(g_pmats(ipm))
        write(cstemp,'(a,i1)') 'pmata_',ipm
        if (.not.nga_create(MT_DBL,2,dims,cstemp(1:7),chunk,
     &     g_pmata(ipm))) call
     &     errquit('dim_addop: nga_create failed '//cstemp(1:7),
     &     0,GA_ERR)
        call ga_zero(g_pmata(ipm))
        write(cstemp,'(a,i1)') 'h1mat_',ipm
        if (.not.ga_create(MT_DBL, nbf, nbf, cstemp(1:7), 0, 0,
     &     g_h1mat(ipm))) call
     &     errquit('dim_addop: ga_create failed '//cstemp(1:7),
     &     0,GA_ERR)
c        if (.not.nga_create(MT_DBL,2,dims,cstemp(1:7),chunk,
c     &     g_h1mat(ipm))) call
c     &     errquit('dim_addop: nga_create failed '//cstemp(1:7),
c     &     0,GA_ERR)
        call ga_zero(g_h1mat(ipm))
      enddo

      dims(1) = 3
      dims(2) = nbf
      dims(3) = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1
      do ipm = 1, 2
        if (.not. nga_create (MT_DBL, 3, dims, 'CPKS dens_r',chunk,
     &     g_dens_r(ipm)))
     &     call errquit('dim_addop: could not allocate g_dens_r',555,
     &     GA_ERR)
        call ga_zero(g_dens_r(ipm))
        if (.not. nga_create (MT_DBL, 3, dims, 'CPKS dens_i',chunk,
     &     g_dens_i(ipm)))
     &     call errquit('dim_addop: could not allocate g_dens_i',555,
     &     GA_ERR)
        call ga_zero(g_dens_i(ipm))
      end do 

      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = 1
      bhi(1) = nbf
      blo(2) = 1
      bhi(2) = nbf
      if (ldebug .and. id .eq. 0) then
        write(luout,*) "End allocation"
      end if
c
      if (.not.rtdb_get(rtdb, 'dimqm:noseed', mt_log, 1, dimqm_noseed))
     $   call errquit('dim_addop: noseed get failed', 1, RTDB_ERR)
c     Pull in residual and check seeding status
      if (.not.lfirst) then
        if (.not.rtdb_get(rtdb, 'lkain:rmax', mt_dbl, 1, rmax))
     $    call errquit('dimqm_addop: rmax get failed', 1, RTDB_ERR)
      
        call dimqm_check_dipoles(rtdb, 1.0d0, rmax)
      end if

c
      do ivec = 1, 3
c
c     ==========================
c     Real Portion of QM Density
c     ==========================
c
        if (ldebug .and. id .eq. 0) then
          write(luout,*) "Start Real " // direction(ivec)
        end if
        do ipm = 1, ncomp
          call ga_zero(g_h1mat(ipm))
          call ga_vec_to_mat(g_h1mat(ipm), 1, nvir, 1, nclosed,
     $         g_x_r(ipm), 1, ivec) ! g_h1mat now holds A-matrix for
          call ga_zero(g_pmats(ipm))
          call ga_zero(g_pmata(ipm))
        end do                  ! ipm = 1,2
        call ga_sync()
c       note: the last argument tells it not to use an occ-occ
c             block to build the density marix.
        call CalcPerturbedTDPmat1_opt
     &     (ncomp, g_pmats, g_pmata, g_h1mat, g_movecs, nbf, nclosed,
     &     nvir, nmo, .false., .false.,
     &     limag, .false.)  ! density matrix -> pmats

c        call ga_zero(g_pmata(1))
c        call ga_zero(g_pmata(2))
        call ga_sync()
        call ga_scale(g_pmats(1),0.25d0)
        call ga_scale(g_pmats(2),0.25d0)
c       
        alo(1) = ivec
        ahi(1) = ivec
        if (.not.limag) then
c       this works for real, symmetric, perturbations
c       calculate P(S) = [P(+) + P(-)]/2
          call nga_add_patch (0.5d0, g_pmats(1), blo, bhi,
     &      0.5d0, g_pmats(2), blo, bhi,
     &      g_dens_r(1), alo, ahi)
c         caluclate P(A) = [-P(+) + P(-)]/2  (wrong results with opposite sign ...)
          call nga_add_patch (-0.5d0, g_pmats(1), blo, bhi,
     &      0.5d0, g_pmats(2), blo, bhi,
     &      g_dens_r(2), alo, ahi)
        else
c
c       this here is for imaginary, antisymmetric, perturbations
c       calculate P(S) = [P(+) - P(-)]/2
          call nga_add_patch (0.5d0, g_pmats(1), blo, bhi,
     &      -0.5d0, g_pmats(2), blo, bhi,
     &      g_dens_r(1), alo, ahi)
c       caluclate P(A) = -[P(+) + P(-)]/2  ! sign needs to be determined
          call nga_add_patch (-0.5d0, g_pmats(1), blo, bhi,
     &      -0.5d0, g_pmats(2), blo, bhi,
     &      g_dens_r(2), alo, ahi)
        end if
        call ga_sync()

        if (ldebug .and. id .eq. 0) then
          write(luout,*) "End Real " // direction(ivec)
        end if
c
c     ===============================
c     Imaginary Portion of QM Density
c     ===============================
c
        if (lifetime) then
          if (ldebug .and. id .eq. 0) then
            write(luout,*) "Start Imag " // direction(ivec)
          end if
          do ipm = 1, ncomp
            call ga_zero(g_h1mat(ipm))
            call ga_vec_to_mat(g_h1mat(ipm), 1, nvir, 1, nclosed,
     $           g_x_i(ipm), 1, ivec) ! g_h1mat now holds A-matrix for
            call ga_zero(g_pmats(ipm))
            call ga_zero(g_pmata(ipm))
          end do                  ! ipm = 1, ncomp
          call ga_sync()
c         note: the last argument tells it not to use an occ-occ
c               block to build the density marix.
          call CalcPerturbedTDPmat1_opt
     &       (ncomp, g_pmats, g_pmata, g_h1mat, g_movecs, nbf,
     &        nclosed, nvir, nmo, .false., .false.,
     &       limag, .false.)  ! density matrix -> pmats
c           write(luout,*) "pmata"
c           call ga_print(g_pmata(1))
c           call ga_print(g_pmata(2))
c zeroing the antisymmetric dens mats is done in get_dens_reorim to force uncoupled CPKS. 
c Not sure if this should be forced for DIM/QM density or not. 
          call ga_zero(g_pmata(1))
          call ga_zero(g_pmata(2))
c looks like g_pmata is not used, so w/e
          call ga_sync()
          call ga_scale(g_pmats(1),0.25d0)
          call ga_scale(g_pmats(2),0.25d0)
c       
          alo(1) = ivec
          ahi(1) = ivec
          if (.not.limag) then
c       calculate P(S) = [P(+) + P(-)]/2
            call nga_add_patch (0.5d0, g_pmats(1), blo, bhi,
     &        0.5d0, g_pmats(2), blo, bhi,
     &        g_dens_i(1), alo, ahi)
c         caluclate P(A) = [-P(+) + P(-)]/2  (wrong results with opposite sign ...)
            call nga_add_patch (-0.5d0, g_pmats(1), blo, bhi,
     &        0.5d0, g_pmats(2), blo, bhi,
     &        g_dens_i(2), alo, ahi)
          else
            call nga_add_patch (0.5d0, g_pmats(1), blo, bhi,
     &        -0.5d0, g_pmats(2), blo, bhi,
     &        g_dens_i(1), alo, ahi)
c       caluclate P(A) = -[P(+) + P(-)]/2  ! sign needs to be determined
c                                          
            call nga_add_patch (-0.5d0, g_pmats(1), blo, bhi,
     &        -0.5d0, g_pmats(2), blo, bhi,
     &        g_dens_i(2), alo, ahi)
          end if ! limag
          if (ldebug .and. id .eq. 0) then
            write(luout,*) "End Imag " //  direction(ivec)
          end if
        end if ! lifetime
        call ga_sync()
      end do ! ivec = 1, 3
c
c   Allocate new arrays
c      if(.not.ma_push_get(mt_dbl,i3*2,'dimrsp:fld',l_fld,k_fld))
c     $  call errquit('malloc dimrsp:fld failed',1,MA_ERR)
c
c
      if(.not.ma_push_get(mt_dbl,i3,'dimrsp:xyz',l_dimxyz,k_dimxyz))
     $  call errquit('malloc dimrsp:xyz failed',1,MA_ERR)
      if(.not.rtdb_get(rtdb,'dimqm:coords',mt_dbl,i3,dbl_mb(k_dimxyz)))
     $  call errquit('dim_addop2: rtdb get failed',1,RTDB_ERR)
      if(.not.ga_create(mt_dbl,i3,1,'g_fld1',0,0,g_fld1))
     $  call errquit('dim_addop2: create g_fld1 failed',0,GA_ERR)
      if(.not.ga_create(mt_dbl,i3,1,'g_fld2',0,0,g_fld2))
     $  call errquit('dim_addop2: create g_fld2 failed',0,GA_ERR)
      if(.not.ga_create(mt_dcpl,i3,1,'g_fldcmplx',0,0,g_fldcmplx))
     $  call errquit('dim_addop2: create g_fldcmplx failed',0,GA_ERR)
      if(.not.ga_create(mt_dcpl,i3,3,'g_fld_tot',0,0,g_fld_tot))
     $  call errquit('dim_addop2: create g_fld_tot failed',0,GA_ERR)
      if(.not.ga_create(mt_dbl,i3,3,'g_mu_r1',0,0,g_mu_r(1)))
     $  call errquit('dim_addop2: create g_mu_r failed',0,GA_ERR)
      if(.not.ga_create(mt_dbl,i3,3,'g_mu_r2',0,0,g_mu_r(2)))
     $  call errquit('dim_addop2: create g_mu_r failed',0,GA_ERR)
      if(.not.ga_create(mt_dbl,i3,3,'g_mu_i1',0,0,g_mu_i(1)))
     $  call errquit('dim_addop2: create g_mu_r failed',0,GA_ERR)
      if(.not.ga_create(mt_dbl,i3,3,'g_mu_i2',0,0,g_mu_i(2)))
     $  call errquit('dim_addop2: create g_mu_r failed',0,GA_ERR)
c
      if(.not.rtdb_get(rtdb,'dimqm:coords', mt_dbl, i3, dimxyz))
     $  call errquit('get dimqm:coords failed', 1, RTDB_ERR)
c
      g_dens_comp_r = ga_create_atom_blocked(geom,basis,
     $                          'real density matrix comp')
      if (lifetime) then
        g_dens_comp_i = ga_create_atom_blocked(geom,basis,
     $                          'imag density matrix comp')
      end if
c
c
c      call dimqm_screening(rtdb, geom, basis, dbl_mb(k_dimxyz),
c     $                       screen)
c      screen = ONE
c
c   =============================
c   Solve for induced dipoles +/-
c   =============================
c
c     Loop over perturbations
      do ipm = 1, 2
        do n = 1, 3
          call ga_zero(g_dens_comp_r)
          if (lifetime) call ga_zero(g_dens_comp_i)
          alo(1) = n
          ahi(1) = n
c
c       Copy current perturbation into g_dens_comp
          call nga_copy_patch('N',g_dens_r(ipm), alo, ahi,
     $                          g_dens_comp_r, blo, bhi)
          if (lifetime) then
            call nga_copy_patch('N',g_dens_i(ipm), alo, ahi,
     $                              g_dens_comp_i, blo, bhi)
          end if
          muind = ZERO
          fld = ZERO
          firsttime = .false.
!          if(.not.rtdb_get(rtdb,
!     $                'dimqm:muind_'//direction(n)//'_r'//pm(ipm),
!     $                              mt_dbl, i3, muold(:,:,1))) then
!            if(id.eq.0 .and. ldebug) 
!     $         write(luout,*) "First cycle, no old dipoles!"
!            muold = ZERO
!            firsttime = .true.
!            dimqm_seeded = .false.
!            call dimqm_put_seeded(rtdb, dimqm_seeded)
!c            xyz_seeded(3*(n-1)+ipm) = .false.
!            if(dimtol0 < 1.0d-4 .and. .not. dimqm_noseed) then
!              dimtolxyz(ipm*3 - 1 + n) = 1.0d-4
!              if(id.eq.0 .and. ldebug) then
!                write(luout,*) "Requested tolerance below 1.0d-4"
!                write(luout,*) "Setting "//direction(n)//pm(ipm)//
!     $                         " dir tolerance to 1.0d-4 to start"
!              end if
!            end if
!          else
!            if(.not.rtdb_get(rtdb,
!     $                'dimqm:muind_'//direction(n)//'_i'//pm(ipm),
!     $                           mt_dbl, i3, muold(:,:,2)))
!     $          call errquit('get dimqm:muold failed',1,RTDB_ERR)
!          end if
c         Set convergence tolerance
c          dimtol = dimtolxyz(ipm*3 - 1 + n)
c          dimqm_seeded = xyz_seeded(ipm*3 - 1 + n)
c          dimtol = 1.0d-7
c          dimqm_noseed = .true.
c          call dfill(i3*2, ZERO, dbl_mb(k_muind), 1)
c          call dfill(i3*2, ZERO, dbl_mb(k_fld), 1)
c
c       Real portion of E-Field
c        write(luout,*) "REAL"
          call dimqm_EqmE_FD(
     $                        rtdb, 
     $                        g_dens_comp_r, 
     $                        geom, 
     $                        basis,
     $                        dbl_mb(k_dimxyz),
     $                        g_fld1) 
c
c       Imaginary portion of E-Field
          call dimqm_EqmE_FD(
     $                        rtdb, 
     $                        g_dens_comp_i, 
     $                        geom, 
     $                        basis,
     $                        dbl_mb(k_dimxyz),
     $                        g_fld2) 
          call dim_convert_d2z(i3,
     $                        g_fld1,
     $                        g_fld2,
     $                        g_fldcmplx)
c
c       Solve for induced dipoles
          call dimqm_iterative_FD(
     $                        rtdb,
     $                        g_fldcmplx,
     $                        dbl_mb(k_dimxyz),
     $                        DIM_nAtoms)
          call ga_copy_patch('n', g_fldcmplx, 1, i3, 1, 1,
     $                            g_fld_tot,  1, i3, n, n)
        end do ! n = 1, 3
        call convert_z2d(g_fld_tot, g_mu_r(ipm), g_mu_i(ipm))
      end do ! ipm = 1, 2
c      if(MAXVAL(dimErr) <= 1.0d-4) then
c        write(luout,*) "Dipole error below 1d-4"
c        write(luout,*) "Shutting down DIM"
c        dimqm_on = .false.
c      end if
c
c   Destroy GAs we don't need anymore
      if (.not. ga_destroy(g_dens_comp_r)) call errquit
     $    ('addop: dens_comp_r GA?',0, GA_ERR)
      if (lifetime) then
      if (.not. ga_destroy(g_dens_comp_i)) call errquit
     $    ('addop: dens_comp_i GA?',0, GA_ERR)
      end if
      do ipm = 1,2
        if (.not. ga_destroy(g_dens_r(ipm))) call errquit
     $     ('addop: dens_r GA?',0, GA_ERR)
        if (.not. ga_destroy(g_dens_i(ipm))) call errquit
     $     ('addop: dens_i GA?',0, GA_ERR)
      end do
c
c
c   ====================================================
c   Solve for DIM potential, both real and imaginary S/A
c   ====================================================
c
      dims(1) = 3
      dims(2) = nbf
      dims(3) = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1
c
c   Real +
      if (.not. nga_create (MT_DBL, 3, dims, 'addop:dim_r+',chunk,
     &                        g_dim_r(1)))
     &   call errquit('addop: could not allocate g_dim_r+',1,GA_ERR)
      call ga_zero(g_dim_r(1)) 
      call dim_grid_FD_gen(
     $                   rtdb,
     $                   g_dim_r(1),
     $                   g_mu_r(1),
     $                   1,
     $                   1)
      call ga_symmetrize(g_dim_r(1))
c
c   Real -
      if (.not. nga_create (MT_DBL, 3, dims, 'addop:dim_r-',chunk,
     &                        g_dim_r(2)))
     &   call errquit('addop: could not allocate g_dim_r-',1,GA_ERR)
      call ga_zero(g_dim_r(2))
      call dim_grid_FD_gen(
     $                   rtdb,
     $                   g_dim_r(2),
     $                   g_mu_r(2),
     $                   1,
     $                   1)
      call ga_antisymmetrize(g_dim_r(2))
      if (lifetime) then
c
c   Imaginary +
      if (.not. nga_create (MT_DBL, 3, dims, 'addop:dim_i+',chunk,
     &                        g_dim_i(1)))
     &   call errquit('addop: could not allocate g_dim_i+',1,GA_ERR)
      call ga_zero(g_dim_i(1))
      call dim_grid_FD_gen(
     $                   rtdb,
     $                   g_dim_i(1),
     $                   g_mu_i(1),
     $                   1,
     $                   1)
      call ga_symmetrize(g_dim_i(1))
c
c   Imaginary -
      if (.not. nga_create (MT_DBL, 3, dims, 'addop:dim_i-',chunk,
     &                        g_dim_i(2)))
     &   call errquit('addop: could not allocate g_dim_i-',1,GA_ERR)
      call ga_zero(g_dim_i(2))
      call dim_grid_FD_gen(
     $                   rtdb,
     $                   g_dim_i(2),
     $                   g_mu_i(2),
     $                   1,
     $                   1)
      call ga_antisymmetrize(g_dim_i(2))
      end if
c
c   ======================================
c   Undo the symmetrization to recover +/-
c   ======================================
c
      do ivec = 1, 3
        alo(1) = ivec
        ahi(1) = ivec
c       ************
c       Real portion
c       ************
        call nga_copy_patch ('N',g_dim_r(1),alo,ahi,g_pmats(1),blo,bhi)
        call nga_copy_patch ('N',g_dim_r(2),alo,ahi,g_pmats(2),blo,bhi)
c
c       it might be necessary to use 0.5 here instead of 1.0
c       (note: that turned out NOT to be the case after some testing)
        pre_factor = 1.0d0
        call ga_sync()
        if (.not.limag) then
c         real perturbation: 
          call nga_add_patch (pre_factor, g_pmats(1), blo, bhi,
     &       pre_factor, g_pmats(2), blo, bhi,
     &       g_dim_r(1), alo, ahi)
          call nga_add_patch (pre_factor, g_pmats(1), blo, bhi,
     &       -pre_factor, g_pmats(2), blo, bhi,
     &       g_dim_r(2), alo, ahi)
        else
c         imaginary perturbation:
          call nga_add_patch (pre_factor, g_pmats(1), blo, bhi,
     &       pre_factor, g_pmats(2), blo, bhi,
     &       g_dim_r(1), alo, ahi)
          call nga_add_patch (-pre_factor, g_pmats(1), blo, bhi,
     &       pre_factor, g_pmats(2), blo, bhi,
     &       g_dim_r(2), alo, ahi)
        end if  ! if .not.limag
        if (lifetime) then
c       *****************
c       Imaginary portion
c       *****************
        call nga_copy_patch ('N',g_dim_i(1),alo,ahi,g_pmats(1),blo,bhi)
        call nga_copy_patch ('N',g_dim_i(2),alo,ahi,g_pmats(2),blo,bhi)
c
c       it might be necessary to use 0.5 here instead of 1.0
c       (note: that turned out NOT to be the case after some testing)
        pre_factor = 1.0d0
        call ga_sync()
        if (.not.limag) then
c         real perturbation: 
          call nga_add_patch (pre_factor, g_pmats(1), blo, bhi,
     &       pre_factor, g_pmats(2), blo, bhi,
     &       g_dim_i(1), alo, ahi)
          call nga_add_patch (pre_factor, g_pmats(1), blo, bhi,
     &       -pre_factor, g_pmats(2), blo, bhi,
     &       g_dim_i(2), alo, ahi)
        else
c         imaginary perturbation:
          call nga_add_patch (pre_factor, g_pmats(1), blo, bhi,
     &       pre_factor, g_pmats(2), blo, bhi,
     &       g_dim_i(1), alo, ahi)
          call nga_add_patch (-pre_factor, g_pmats(1), blo, bhi,
     &       pre_factor, g_pmats(2), blo, bhi,
     &       g_dim_i(2), alo, ahi)
        end if  ! if .not.limag
        end if ! lifetime
      enddo                     ! ivec = 1,nvec
100   continue
c
c   ====================================
c   Add DIM potential to the Fock matrix
c   ====================================
c
      g_dcv = ga_create_atom_blocked(geom, basis, 'rohf_h2e3: dcv')
      xoff = 1
      voff = nclosed + nopen + 1
      xend = nvir*nclosed
      do ivec = 1, 3 ! Loop over perturbations
        alo(1) = ivec
        ahi(1) = ivec
        do ipm = 1, ncomp! Loop over +/-
c         We only add the + direction of the DIM potential to both +/- of the Fock matrix
c   Real Portion
          call nga_copy_patch('N',g_dim_r(ipm),alo,ahi,g_dcv,blo,bhi)
          call ga_scale(g_dcv, two)
          call ga_matmul_patch('n', 'n', two, zero,
     $                           g_dcv,   1, nbf, 1, nbf,
     $                           g_movecs, 1, nbf, 1, nclosed,
     $                           g_tmp1,  1, nbf, 1, nclosed)
          call ga_sync()
          call ga_matmul_patch('t', 'n', one, zero,
     $                           g_movecs, voff, nmo, 1, nbf,
     $                           g_tmp1, 1, nbf,  1, nclosed,
     $                           g_tmp2, 1, nvir, 1, nclosed)
          call ga_sync()
          call ga_mat_to_vec(g_tmp2, 1, nvir, 1, nclosed, g_Ax_r(ipm),
     $                         xoff, ivec, four, '+')
c
c   Imaginary Portion
          if (lifetime) then
          call nga_copy_patch('N',g_dim_i(ipm),alo,ahi,g_dcv,blo,bhi)
          call ga_scale(g_dcv, two)
          call ga_matmul_patch('n', 'n', two, zero,
     $                           g_dcv,   1, nbf, 1, nbf,
     $                           g_movecs, 1, nbf, 1, nclosed,
     $                           g_tmp1,  1, nbf, 1, nclosed)
          call ga_sync()
          call ga_matmul_patch('t', 'n', one, zero,
     $                           g_movecs, voff, nmo, 1, nbf,
     $                           g_tmp1, 1, nbf,  1, nclosed,
     $                           g_tmp2, 1, nvir, 1, nclosed)
          call ga_sync()
c
c   ******NOTE JEM******
c   We can remove the - sign if we apply the sign change discussed in aoresponse and rohf_hessv_xx3
c   ********************
c
          call ga_mat_to_vec(g_tmp2, 1, nvir, 1, nclosed, g_Ax_i(ipm),
     $                         xoff, ivec, four, '+')
          end if !lifetime
        end do !ipm = 1, 2
      end do !ivec = 1, 3
c   ========
c   Clean up
c   =======
      lfirst = .false.
      call dimqm_setlfirst(rtdb, lfirst)
      stat = .true.
      stat = stat.and.ma_chop_stack(l_dimxyz)
      stat = stat.and.ga_destroy(g_fld1)
      stat = stat.and.ga_destroy(g_fld2)
      stat = stat.and.ga_destroy(g_fldcmplx)
      stat = stat.and.ga_destroy(g_fld_tot)
      stat = stat.and.ga_destroy(g_mu_r(1))
      stat = stat.and.ga_destroy(g_mu_r(2))
      stat = stat.and.ga_destroy(g_mu_i(1))
      stat = stat.and.ga_destroy(g_mu_i(2))
      if (.not. stat) 
     $    call errquit('dimqm_addop2 cleanup failed',0,GA_ERR)

      do ipm = 1,2
        if (.not. ga_destroy(g_pmats(ipm))) call errquit
     $     ('addop: pmats GA?', 0, GA_ERR)
        if (.not. ga_destroy(g_pmata(ipm))) call errquit
     $     ('addop: pmata GA?', 0, GA_ERR)
        if (.not. ga_destroy(g_h1mat(ipm))) call errquit
     $     ('addop: h1mat GA?', 0, GA_ERR)
        if (.not. ga_destroy(g_dim_r(ipm))) call errquit
     $     ('addop: dim_r GA?', 0, GA_ERR)
        if (lifetime) then
        if (.not. ga_destroy(g_dim_i(ipm))) call errquit
     $     ('addop: dim_i GA?', 0, GA_ERR)
        end if
      enddo                     ! ipm = 1,2
c
      if (.not. ga_destroy(g_tmp1)) call errquit
     $   ('addop: tmp1 GA?', 0, GA_ERR)
      if (.not. ga_destroy(g_tmp2)) call errquit
     $   ('addop: tmp2 GA?', 0, GA_ERR)
      if (.not. ga_destroy(g_dcv)) call errquit
     $   ('addop: dcv GA?',0, GA_ERR)
c
      end subroutine dimqm_addop2
