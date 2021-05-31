      subroutine dimqm_gs(rtdb, basis, geom, g_dens, lastcy)
c     Author: Jeff Becca
c     called from: nwdft/scf_dft/dft_scf.F
c     subroutines called:  dimqm_eqmn_gs, dimqm_eqme_gs, dimqm_f2d2
c ----------------------------------------------------------------------
c  This is the main driver for DIM/QM. It is called every ground state
c  SCF cycle. currently modeled off of dimqm_main
c ----------------------------------------------------------------------
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "dimqm.fh"
c------------------------- Input Variables -----------------------------
      integer rtdb      ! RTDB handle
      integer basis     ! Basis set handle
      integer geom      ! Geom handle
      integer g_dens    ! QM density handle (SHOULD BE TOTAL DENS HERE)
      logical lastcy    ! whether or not this is the last DIM cycle
c------------------------ Local Variables ------------------------------
      logical ldebug
      logical stat
      logical lfirst
      integer id, i, i3
      integer l_dimxyz, k_dimxyz
      integer l_eqmn, k_eqmn
      integer l_eqme, k_eqme
      integer l_etot, k_etot
C      integer l_muold, k_muold
      double precision muold(DIM_nOrder)
      integer l_muind, k_muind
      double precision edimqm, dimqm_energy2, edimqmN, edimqmE
      double precision dimtol
      double precision maxdiperr
      character*50 d
      character*50 dd
c----------- GA's -----------
      integer g_eqme, g_eqmn, g_etot
c----------------------------
c| Initialize some settings |
c----------------------------
      id = ga_nodeid()
      stat = rtdb_parallel(.true.)
      call dimqm_getlfirst(rtdb, lfirst)

      if (lfirst) call ga_zero(g_DIMmuold)
c ----------------- Initialize GA's used in DIM cycle ------------------
      if(.not.ga_create(mt_dbl,DIM_nOrder,1,'g_eqme',0,0,g_eqme))
     $   call errquit('dimqm_gs: create g_eqme failed',0,GA_ERR)
      call ga_zero(g_eqme)
      if(.not.ga_create(mt_dbl,DIM_nOrder,1,'g_eqmn',0,0,g_eqmn))
     $   call errquit('dimqm_gs: create g_eqmn failed',0,GA_ERR)
      call ga_zero(g_eqmn)
      if(.not.ga_create(mt_dbl,DIM_nOrder,1,'g_etot',0,0,g_etot))
     $   call errquit('dimqm_gs: create g_etot failed',0,GA_ERR)
      call ga_zero(g_etot)
c  allocate a bunch of persistant arrays. Kill at the end of main
      if(.not.ma_push_get(mt_dbl,3*DIM_nAtoms,'dimqm:dimxyz',l_dimxyz,
     $   k_dimxyz)) call errquit('dimqm_gs: malloc k_dimxyz failed',
     $   1,MA_ERR)
      if(.not.ma_alloc_get(mt_dbl,3*DIM_nAtoms,'dimqm:eqme',l_eqme,
     $   k_eqme)) call errquit('dimqm_gs: malloc k_eqme failed',
     $   1,MA_ERR)
      if(.not.ma_alloc_get(mt_dbl,3*DIM_nAtoms,'dimqm:muind',l_muind,
     $   k_muind)) call errquit('dimqm_gs: malloc k_muind failed',
     $   1,MA_ERR)
c  get and allocate coords
      if(.not.rtdb_get(rtdb,'dimqm:coords',mt_dbl,
     $      3*DIM_nAtoms,dbl_mb(k_dimxyz))) call errquit
     $      ('dimmain2: xyz rtdb get failed',1,RTDB_ERR)
c  Call EqmN2 to get the field contributions due to QM nuclei
      call dimqm_eqmn_gs(rtdb, geom, dbl_mb(k_dimxyz), g_eqmn)
c  Zero EqmE
      call dfill(3*DIM_nAtoms, 0.0d0, dbl_mb(k_eqme), 1)
c  Call EqmE to get field contributions due to QM e-
      call dimqm_eqme_gs(rtdb, g_dens, geom, basis, 
     $         dbl_mb(k_eqme), dbl_mb(k_dimxyz), DIM_nAtoms, g_eqme)
c  put last cycles dipoles into g_DIMmuold and zero g_DIMmuind for this cycle
      call ga_copy(g_DIMmuind, g_DIMmuold)
      call ga_zero(g_DIMmuind)
c  Combine field sources, store in g_DIMmuind
      !call ga_add(1.0d0, g_eqme, 1.0d0, g_eqmn, g_etot)
      call ga_add(1.0d0, g_eqme, 1.0d0, g_eqmn, g_DIMmuind)
c  ------------------------- Things we need from calc ------------------
c     Name           reason
c     ----           ------   
c  g_muind           ind. dipoles due to total fields
c  g_eqme            int. energy due to e-
c  g_eqmn            int. energy due to nu+
c  g_muold           last cycle dipoles to see if dipoles converged
c  g_etot            total E-fields for calculating ind. dip.

c  ** g_muold does not need to be passed to anything, we either have 
c     them or dont **
c
c  Only thing f2d2 should need is array for ind. dip. (g_muind),
c  total Efield (g_etot), and DIM coordinates (k_dimxyz)
c  
c  ---------------------------------------------------------------------
      call dimqm_iterative2(rtdb, g_DIMmuind, 
     $                        dbl_mb(k_dimxyz), DIM_nAtoms)
ccc Calculate interaction energy and put it in rtdb
      edimqmN = -0.5d0 * ga_ddot(g_DIMmuind, g_eqmn)
      edimqmE = -0.5d0 * ga_ddot(g_DIMmuind, g_eqme)
      edimqm = edimqmN + edimqmE

      if(.not.rtdb_put(rtdb, 'dimqm:edimqm',mt_dbl,1,edimqm))
     $   call errquit('main2: edimqm put failed',1,RTDB_ERR)
c  TODO: store ind. dip. in rtdb for now, so they can be accessed
c        when creating potential
c        ** NOTE ** it should be looked into for having the potential 
c        created and stored to rtdb instead. should be much smaller
c        than the dipoles since it only needs to apply to QM grid   
c        ** NEWER NOTE ** screw this writing to disk nonsense. Lets
c        keep it all in memory. 
      call ga_init_fence()
      call ga_get(g_DIMmuind, 1, DIM_nOrder, 1, 1, dbl_mb(k_muind), 1)
      call ga_fence()
c      call ga_print(g_etot)
c      if (id.eq.0) then
c         do i = 0, nOrder-1
c            write(luout,*)'k_muind',i+1,' ',dbl_mb(k_muind+i)
c         enddo
c      endif
      if(.not.rtdb_put(rtdb, 'dimqm:muind',mt_dbl,DIM_nOrder,
     $   dbl_mb(k_muind))) call errquit('dim_main2: rtdb put muind
     $   failed',1,RTDB_ERR)
c  TODO: Check that dipoles are converged
      call ga_init_fence()
      call ga_get(g_DIMmuold, 1, DIM_nOrder, 1, 1, muold, 1)
      call ga_fence()
      call dimqm_dip_err(dbl_mb(k_muind), muold, maxdiperr)

c Place this cycle's dipoles into g_muold
      call ga_copy(g_DIMmuind, g_DIMmuold)

cc TODO  gather some nice info for printing out the induced dipoles TODO
cc         
cc      maxerr = MAXVAL(ABS(muold-muind))
c      call dimqm_dip_err(dbl_mb(k_muold), dbl_mb(k_muind),nDIM,maxerr)
c      dimtol = 1.0d-5
c      call dimqm_check_dipoles(rtdb, dimtol, maxerr)
cc
ccc     make this cycle's dipoles the old dipoles for next cycle
ccc      muold = 0.0d0
c      call dfill(3*nDIM, 0.0d0, dbl_mb(k_muold), 1)
ccc      muold = muind
cc      do i = 1, nDIM
cc         dbl_mb(k_muold+(i-1)) = dbl_mb(k_muind+(i-1))
cc      enddo
      call dimqm_setlfirst(rtdb,.false.)
      stat = .true.
      stat = stat.and.ma_free_heap(l_eqme)
      stat = stat.and.ma_free_heap(l_muind)
      stat = stat.and.ma_chop_stack(l_dimxyz)
      stat = stat.and.ga_destroy(g_eqme)
      stat = stat.and.ga_destroy(g_eqmn)
      stat = stat.and.ga_destroy(g_etot)
      if (.not.stat) 
     $   call errquit('dimqm_gs:free stack failed',1,MA_ERR)


      dd =
     $ '=================================================='
      d =
     $ '--------------------------------------------------'
      if(lastcy) then
         if(id.eq.0) then
c  TODO
            write(luout,*)
     $'---------------------------- DIM/QM RESULTS -------------------'
            write(luout,*)''

            write(LuOut,*) 'DIM/QM Energy'
            write(LuOut,*) d
            write(LuOut,111) edimqmE
            write(LuOut,222) edimqmN
c      write(LuOut,444) efld_int
            write(LuOut,333) edimqm
         endif

         call dimqm_printAtoms(rtdb, dbl_mb(k_dimxyz),
     $                           DIM_nAtoms,DIM_nTypes)

        call dimqm_printDipole(rtdb, g_DIMmuind, 1, .true., DIM_nAtoms)
         call dimqm_destroy(rtdb)
      endif

 111  format(' Polarization (Dipole, el.) = ',f20.10)
 222  format(' Polarization (Dipole, nuc) = ',f20.10)
 333  format('              Total Energy  = ',f20.10)
 444  format(' Interaction with Efield    = ',f20.10)
 555  format('                      Total = ',f20.10)


      end subroutine dimqm_gs

      subroutine dimqm_dip_err(new, old, maxerr)
c-----------------------------------------------------------------------
c  Author: Jeff Becca
c  Calls: Nothing
c  Called from: dimqm_gs
c  Purpose: get max err from allocated arrays
      implicit none
#include "global.fh"
#include "stdio.fh"
#include "dimqm.fh"
c------------------------------ Input variables ------------------------
      integer a, b
      double precision maxerr, new(DIM_nOrder), old(DIM_nOrder)

      maxerr = MAXVAL(ABS(old - new))
      write(luout,*)'max dip error: ', maxerr
      end subroutine dimqm_dip_err

      function dimqm_energy2(muind, eqme, eqmn, DIM_nAtoms)
c-----------------------------------------------------------------------
c  Author: Jeff Becca
c  Calls: Nothing
c  Called from: dimqm_gs
c  Purpose: Calculate DIM/QM interaction energy
c-----------------------------------------------------------------------
      implicit none
c-------------------------- Input variables ----------------------------
      double precision muind(3*DIM_nAtoms)
      double precision  eqme(3*DIM_nAtoms)
      double precision  eqmn(3*DIM_nAtoms)
      double precision dimqm_energy2, tmp1, tmp2
      integer DIM_nAtoms, i

c      tmp1 = -0.5d0 * SUM(muind * eqme)
c      tmp2 = -0.5d0 * SUM(muind * eqmn)
c      dimqm_energy2 = tmp1 + tmp2
      tmp1 = 0.0d0
      tmp2 = 0.0d0
      do i = 1, 3*DIM_nAtoms
         tmp1 = tmp1 + (muind(i)*eqme(i))
         tmp2 = tmp2 + (muind(i)*eqmn(i))
      enddo
      dimqm_energy2 = -0.5d0*(tmp1 + tmp2)
      return
      end function dimqm_energy2

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
