      subroutine dimqm_main(rtdb, basis, geom, g_dens, lastcy, nDIM)
c
c     Author: Justin Moore
c             Jeff Becca
c
c     Called from: dft_scf.F, rohf_fock.F 
c
c     Subroutines called: dimqm_prep, dimqm_EqmE, dimqm_EqmN, dimqm_f2d
c                         dimqm_seed_init_tolerance, dimqm_check_dipoles,
c                         dimqm_energy, dimqm_printAtoms, dimqm_printDipole,
c                         dimqm_printAtomicDipoles, dimqm_tidy 
c
c     Main driver for the DIM/QM module.  Called every SCF cycle.
c    
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
c   
c     Input Variables
      integer rtdb      ! RTDB handle
      integer basis     ! Basis set handle
      integer geom      ! Geometry handle
      integer g_dens    ! QM density GA handle (we only work with the total density
      logical lastcy    ! Last SCF iteration flag
      integer nDIM
c
c     Local variables
      integer i3, m, i                  ! Looping and convenience variables
      double precision EqmE(3, nDIM)    ! Electric field from QM electrons felt at each DIM atom
      double precision EqmN(3, nDIM)    ! Ele. field from QM Nuc. felt at each DIM atom
      double precision Etot(3, nDIM)    ! Total QM electric field (electrons + nuclei) felt at each DIM atom
      integer l_fin,  k_fin             ! Finite field
      double precision dimxyz(3, nDIM)  ! DIM coordinates
      double precision muind(3, nDIM)   ! DIM induced dipoles
      double precision muold(3, nDIM)
      double precision max_err, calcErr ! Error in dipoles
      logical stat                      ! Dummy for function returns
      integer id                        ! Processor ID
      double precision dime             ! DIM/QM interaction energy function
      double precision edimqm
c      integer a
      double precision dx, dy, dz       ! Total induced dipoles moment for the DIM system
      double precision efld_int
      logical lfirst
      logical lclfld
      logical lefield
      logical ldebug
      logical latmdip
      integer nDIMTypes
      double precision dimtol
      double precision finiteFieldX, finiteFieldY, finiteFieldZ
c
c   ===========
c   Initialize
c   =========== 
      i3 = 3 * nDIM
      id = ga_nodeid()
      call dimqm_get_debug(rtdb, ldebug)
      call dimqm_getlefield(rtdb, lefield)
      call dimqm_getlclfld(rtdb, lclfld)
      call dimqm_getlatmdip(rtdb, latmdip)
      write(luout,*)'ldebug',ldebug
      if(id.eq.0.and.ldebug) then
        write(luout,*) "Entering DIM/QM Main"
      end if
c
c   ================
c   Allocate arrays
c   ================
c
      EqmE  = ZERO
      muind = ZERO
      Etot  = ZERO
      if(id.eq.0.and.ldebug) then
        write(luout,*) "Allocated local arrays"
      end if
c
c      Read in coords from RTDB
      if(.not. rtdb_get(rtdb,'dimqm:coords',mt_dbl,i3,dimxyz))
     $   call errquit('dimqm:xyz rtdb get failed',1,RTDB_ERR)
      if(id.eq.0.and.ldebug) then
        write(luout,*) "Pulled coords"
c        write(*,*) "GEOM:", geom
      end if
cc     Read in whether or not this is first cycle
      call dimqm_getlfirst(rtdb, lfirst)
cc
cc   ===========
cc   First Cycle
cc   ===========
cc
      if(lfirst) then
c
c        Allocate memory for QM nuclear contribution on the heap (PERSISTANT)
c        if(.not.ma_alloc_get(mt_dbl,i3,'dimqm eqmn',l_eqmn,k_eqmn))
c     $    call errquit('dimqm_main malloc k_eqmn failed', 1, MA_ERR)
c
c        Allocate memory for last cycle's dipoles on the heap (PERSISTANT)
c        if(.not.ma_alloc_get(mt_dbl,i3,'dimqm muold',l_muold,k_muold))
c     $    call errquit('dimqm_main malloc k_muold failed',1,MA_ERR)
c
        if(id.eq.0.and.ldebug) then
          write(luout,*) "Allocated heap arrays"
        end if
c        Zero both arrays
c        call dfill(i3, ZERO, dbl_mb(k_muold), 1)
        call dfill (i3, ZERO, muold, 1)
c        muold = ZERO
c        call dfill(i3, ZERO, dbl_mb(k_eqmn), 1)
        call dfill(i3, ZERO, EqmN, 1)
c        Build interaction matrix
c        call dimqm_prep(rtdb, -ONE)
c        Generate electric field due to QM nuclei.  This doesn't change so we only do it once.
        write(luout,*)'arrays zerod'
        write(luout,*)'nDIM',nDIM
        call dimqm_EqmN(rtdb, geom, dimxyz,eqmn,nDIM)
        write(luout,*)'back from eqmn'
        if(.not.rtdb_put(rtdb,'dimqm:EqmN',mt_dbl,i3, eqmn))
     $   call errquit('dimqm_main eqmn put failed', 1, RTDB_ERR)
        write(luout,*)'stored eqmn'
      
        call dimqm_seed_init_tolerance(rtdb,1.0d-4)
        if(id.eq.0.and.ldebug) then
          write(luout,*) "DIM/QM First Cycle Prep Done"
        end if
      end if
c
c      Generate E-field from QM electrons
      write(luout,*)'calling EqmE'
      call dimqm_EqmE(rtdb, g_dens, geom, basis, EqmE,
     $                dimxyz, nDIM)
      write(luout,*)'back from EqmE'
c
c      Add finite field if present
      if(lefield .and. lclfld) then
c        if(id.eq.0 .and. .not. lastcy) then 
c          write(luout,*) "Adding finite field to Eqm"
c          write(luout,*) finiteFieldX, finiteFieldY, finiteFieldZ
c        end if
        if(.not.rtdb_get(rtdb,'dimqm:finiteFieldX',mt_dbl,
     $   1,finiteFieldX)) call errquit('get fieldx failed',1,RTDB_ERR)
        if(.not.rtdb_get(rtdb,'dimqm:finiteFieldY',mt_dbl,
     $   1,finiteFieldY)) call errquit('get fieldy failed',1,RTDB_ERR)
        if(.not.rtdb_get(rtdb,'dimqm:finiteFieldZ',mt_dbl,
     $   1,finiteFieldZ)) call errquit('get fieldz failed',1,RTDB_ERR)
        if(.not.ma_push_get(mt_dbl,i3,'dimqm finite',l_fin,k_fin))
     $   call errquit('dimqm_main malloc k_finite failed',1,MA_ERR)
        call dfill(nDIM, finiteFieldX, dbl_mb(k_fin)  , 3)
        call dfill(nDIM, finiteFieldY, dbl_mb(k_fin+1), 3)
        call dfill(nDIM, finiteFieldZ, dbl_mb(k_fin+2), 3)
c
c        Add finite field for finite-difference local field calculations
c  jbecca TODO: make a separate variable for the finitefield so 
c              the energy can be separated out for static local field
        call daxpy(i3, ONE, dbl_mb(k_fin), 1, EqmE, 1)
        if(.not. ma_pop_stack(l_fin)) call
     $       errquit('ma_pop failed l_fin',1,MA_ERR)
      end if
c
c      Add electronic and nuclear efields into etot
      Etot = EqmE
c      call daxpy(i3, ONE, dbl_mb(k_eqmn), 1, Etot, 1)
      if (.not. lfirst) then
         if(.not.rtdb_get(rtdb, 'dimqm:EqmN',mt_dbl,i3, eqmn))
     $      call errquit('dimqm_main EqmN not read',1,RTDB_ERR)
      endif
      
c      call daxpy(i3, ONE, dbl_mb(k_eqmn), 1, Etot, 1)
      call daxpy(i3, ONE, eqmn, 1, Etot, 1)
c      
c   =======================================================
c   Generate induced dipoles (and charges) on the DIM atoms
c   =======================================================
c
c     since not a common anymore, we need to read old dipoles and place
c     them into muind or it will be very fishy
      if (.not. lfirst) then
         if(.not.rtdb_get(rtdb, 'dimqm:muind',mt_dbl,i3, muold))
     $      call errquit('dimqm_main old dipoles not read',1,RTDB_ERR)
      endif

c      if (ldebug) write(luout,*)'muoldbefore f2d',muold
      if (lefield .and. lclfld) then
         call dimqm_f2d(rtdb, Etot, muind, muold, dimxyz, 1,
     $               'g', nDIM)
      else
         call dimqm_f2d(rtdb, Etot, muind, muold, dimxyz, 1,
     $               'g', nDIM)
      endif
c
      if (ldebug) write(luout,*)'i got here in main'
c      Write induced dipoles and charges to disk
      if(.not.rtdb_put(rtdb,'dimqm:muind',mt_dbl,i3, muind))
     $  call errquit('dimqm_main muind put failed', 1, RTDB_ERR)
c
      if (ldebug) write(luout,*)'i got here in main2'
c      Calculate DIM/QM interaction energy
c      This is needed in the RTDBdd
c      edimqm = dime(EqmE, dbl_mb(k_eqmn), muind, nDIM)
      edimqm = dime(EqmE, eqmn, muind, nDIM)
      if (ldebug) write(luout,*)'edimqm',edimqm

c  jbecca START: need to add in the interaction energy with the efield
c              if the efield is present
c  NOTE: This is not needed for now, since the finitefield is added to
c        EqmE, the interaction energy will also be calculated as part
c        of -HALF * SUM(EqmE * muind)
c
c      if (lefield .and. lastcy) then
c         write(luout,*)'adding in DIM/QM interaction energy with efield'
c         efld_int    =  ZERO
c         efld_int    =  SUM(muind(1,:) * finiteFieldX)
c     $               +  SUM(muind(2,:) * finiteFieldY)
c     $               +  SUM(muind(3,:) * finiteFieldZ)
c         efld_int    =  -HALF * efld_int
c         write(luout,*)'edimqm before efld_int',edimqm
c         write(luout,*)'efld_int',efld_int
c         edimqm      =  edimqm + efld_int
c         write(luout,*)'edimqm after efld_int',edimqm
c      endif       !lefield
c  jbecca END

      if (ldebug) write(luout,*)'i got here in main3'
      if (.not. rtdb_put(rtdb, 'dimqm:edimqm', MT_DBL, 1, edimqm))
     $  call errquit('dimqm_main edimqm put failed', 1, RTDB_ERR)
c
c
      if (ldebug) write(luout,*)'i got here in main4'
      dx = SUM(muind(1,:))
      dy = SUM(muind(2,:))
      dz = SUM(muind(3,:))
      if (ldebug) write(luout,*)'i got here in main5'
      if(id.eq.0 .and. ldebug .and. .not. lastcy) then
        write(luout,*) "Induced dipoles:"
        write(luout,*) dx, dy, dz
        write(luout,*) ''
      end if
      if (ldebug) write(luout,*)'i got here in main6'
c
c      Determine error in dipoles
c      max_err = calcErr(i3, dbl_mb(k_muold), muind)
      max_err = MAXVAL(ABS(muold-muind))
      if (ldebug) write(luout,*)'i got here in main7'
      if(id.eq.0 .and. ldebug .and. .not. lastcy) then
        write(LuOut,*) "Max error in dipoles:", max_err
        write(LuOut,*) ''
      end if
c      if (ldebug) write(luout,*)'muold',muold
c      if (ldebug) write(luout,*)'muind',muind
      if (ldebug) write(luout,*)'i got here in main8'
c     call dimqm_check_dipoles(2.0d-1, max_err)
      if (.not. rtdb_get(rtdb, 'dimqm:dimtol',mt_dbl,1,dimtol))
     $   call errquit('dimqm_main:get dimtol failed',1,RTDB_ERR)
      call dimqm_check_dipoles(rtdb, dimtol, max_err)
c
c      Copy new dipoles into muold
c      call dcopy(i3, muind, 1, dbl_mb(k_muold), 1)
c      call dcopy(i3, muind, 1, muold, 1)
      muold = ZERO
      muold = muind
c      if (ldebug) write(luout,*)'muold2',muold

      if(id.eq.0.and.ldebug) write(LuOut,*)
     $  "Saved dipoles as old dipoles"
c
c      No longer the first cycle!
      lfirst = .false.
      if (ldebug) write(luout,*)'i got here in main9'
      call dimqm_setlfirst(rtdb, lfirst)
c      if (.not. rtdb_put(rtdb, 'dimqm:lfirst', mt_log, 1, lfirst))
c     $   call errquit('dimqm_main: lfirst put failed', 1, RTDB_ERR)
      if (ldebug) write(luout,*)'i got here in main10'
c
c      Last cycle only
      if(lastcy) then
        if(id.eq.0) then
          write(luout,*)
     $    "                                    DIM/QM Results"
          write(luout,*)
     $    "                                    --------------"
          write(luout,*) ""
          write(luout,*) ""

c         write(luout,*)'ndimtypes1', nDIMTypes
         if(.not.rtdb_get(rtdb,'dimqm:ndimtypes',mt_int,1,nDIMTypes))
     $   call errquit('dimqm_main:get dimtypes failed',1,RTDB_ERR)
         write(luout,*)'ndimtypes2', nDIMTypes
c        Calculate energy terms
c        call dimqm_energy(nDIM, muind, EqmE, dbl_mb(k_eqmn))
        call dimqm_energy(rtdb, nDIM, muind, EqmE, eqmn)
c        Print coordinates
        call dimqm_printAtoms(rtdb, dimxyz,nDIM, nDIMTypes)
c        Print total dipole moment
        call dimqm_printDipole(rtdb, muind, 1, .true., nDIM)
c        Print atomic dipoles (and charges)
        if (latmdip) then
            call dimqm_printAtomicDipoles(rtdb, .false., 
     $                                .false., '',nDIM, nDIMTypes)
        endif
c        Destroy persistant arrays
c        if(.not.ma_free_heap(l_muold))
c     $    call errquit("dimqm_main: can't free l_muold", MA_ERR)
c        if(.not.ma_free_heap(l_eqmn))
c     $    call errquit("dimqm_main: can't free l_eqmn", MA_ERR)
c        call dimqm_tidy()
        end if
      end if
c
      if(id.eq.0.and.ldebug) write(LuOut,*) "End dimqm_main"
c      call util_flush(LuOut)
c      call ga_sync()
      end subroutine dimqm_main

      function calcErr(n, old, new)
      implicit none
#include "stdio.fh"
      integer n
      double precision old(n)
      double precision new(n)
      double precision calcErr
      double precision temp(n)
c
      temp(:) = new(:) - old(:)
      calcErr = MAXVAL(temp)
      if (ABS(MINVAL(temp)) > calcErr) then
         calcErr = ABS(MINVAL(temp))
      end if
      return
c
      end function calcErr

      function dime(eqme, eqmn, muind, nAtoms)
      implicit none
#include "dimqm_constants.fh"
      integer nAtoms
      double precision eqme(3,nAtoms)
      double precision eqmn(3,nAtoms)
      double precision muind(3,nAtoms)
      double precision dime
      
      dime = -HALF * SUM(muind * eqme)
      dime = dime +(-HALF)* SUM(muind * eqmn)
      return
      end function dime

      subroutine print_dim(array, nel, label)
      implicit none
#include "stdio.fh"
      integer nel
      double precision array(3,nel)
      character*4 label
      write(luout,*) label
      write(luout,*) array
      call util_flush(LuOut)
      end subroutine print_dim
