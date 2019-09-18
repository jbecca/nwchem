      subroutine dimqm_main2(rtdb, basis, geom, g_dens, lastcy,
     $                           nDIM, g_muold)
c     Author: Jeff Becca
c     called from: nwdft/scf_dft/dft_scf.F
c     subroutines called:  dimqm_EqmN2, dimqm_EqmE2, dimqm_f2d2
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
c------------------------- Input Variables -----------------------------
      integer rtdb      ! RTDB handle
      integer basis     ! Basis set handle
      integer geom      ! Geom handle
      integer g_dens    ! QM density handle (SHOULD BE TOTAL DENS HERE)
      logical lastcy    ! whether or not this is the last DIM cycle
      integer nDIM      ! number of DIM atoms in calculation
      integer g_muold
c------------------------ Local Variables ------------------------------
      logical ldebug
      logical stat
      logical lfirst
      integer id, i, i3, nOrder
      integer l_dimxyz, k_dimxyz
      integer l_eqmn, k_eqmn
      integer l_eqme, k_eqme
      integer l_etot, k_etot
      integer l_muold, k_muold
      integer l_muind, k_muind
      double precision edimqm, dimqm_energy2
      double precision maxerr
      double precision dimtol
c----------- GA's -----------
      integer g_muind, g_eqme, g_eqmn, g_etot
c----------------------------
c| Initialize some settings |
c----------------------------
c TODO: use global arrays as much as possible. Leave xyz as memalloc
c     since dimxyz is not changed, it should be fine
c  since we want the ability to read old dipoles every cycle, I think
c  muold needs to be a GA that is passed to this subroutine. not sure
c  of a way to write GA to disk and read easily
      id = ga_nodeid()
      stat = rtdb_parallel(.true.)
      call dimqm_get_debug(rtdb, ldebug)
      if (ldebug.and.(id.eq.0)) write(luout,*)'starting dimqm_main2'
      call dimqm_getlfirst(rtdb, lfirst)
      nOrder = 3*nDIM
      if (lfirst) call ga_zero(g_muold)
c ----------------- Initialize GA's used in DIM cycle ------------------
c      if(.not.ga_create(mt_dbl,3,nDIM,'g_muind',0,0,g_muind))
c     $   call errquit('dimqm_main2: create g_muind failed',0,GA_ERR)
c      call ga_zero(g_muind)
      if(.not.ga_create(mt_dbl,nOrder,1,'g_eqme',0,0,g_eqme))
     $   call errquit('dimqm_main2: create g_eqme failed',0,GA_ERR)
      call ga_zero(g_eqme)
      if(.not.ga_create(mt_dbl,nOrder,1,'g_eqmn',0,0,g_eqmn))
     $   call errquit('dimqm_main2: create g_eqmn failed',0,GA_ERR)
      call ga_zero(g_eqmn)
      if(.not.ga_create(mt_dbl,nOrder,1,'g_etot',0,0,g_etot))
     $   call errquit('dimqm_main2: create g_etot failed',0,GA_ERR)
      call ga_zero(g_etot)
c  allocate a bunch of persistant arrays. Kill at the end of main
      if(.not.ma_push_get(mt_dbl,3*nDIM,'dimqm:dimxyz',l_dimxyz,
     $   k_dimxyz)) call errquit('dimqm_main2: malloc k_dimxyz failed',
     $   1,MA_ERR)
c      if(.not.ma_alloc_get(mt_dbl,3*nDIM,'dimqm:eqmn',l_eqmn,
c     $   k_eqmn)) call errquit('dimqm_main2: malloc k_eqmn failed',
c     $   1,MA_ERR)
      if(.not.ma_alloc_get(mt_dbl,3*nDIM,'dimqm:eqme',l_eqme,
     $   k_eqme)) call errquit('dimqm_main2: malloc k_eqme failed',
     $   1,MA_ERR)
      if (ldebug.and.(id.eq.0)) write(luout,*)'got here 2 in main'
c      if(.not.ma_alloc_get(mt_dbl,3*nDIM,'dimqm:etot',l_etot,
c     $   k_etot)) call errquit('dimqm_main2: malloc k_etot failed',
c     $   1,MA_ERR)
c      if(.not.ma_alloc_get(mt_dbl,3*nDIM,'dimqm:muold',l_muold,
c     $   k_muold)) call errquit('dimqm_main2: malloc k_muold failed',
c     $   1,MA_ERR)
      if(.not.ma_alloc_get(mt_dbl,3*nDIM,'dimqm:muind',l_muind,
     $   k_muind)) call errquit('dimqm_main2: malloc k_muind failed',
     $   1,MA_ERR)
c  get and allocate coords
         if(.not.rtdb_get(rtdb,'dimpar:coords',mt_dbl,
     $      3*nDIM,dbl_mb(k_dimxyz))) call errquit
     $      ('dimmain2: xyz rtdb get failed',1,RTDB_ERR)
c  Call EqmN2 to get the field contributions due to QM nuclei
      call dimqm_EqmN2(rtdb, geom, dbl_mb(k_dimxyz), g_eqmn, nDIM)
c      call ga_print(g_eqmn)
c  Zero EqmE
      call dfill(3*nDIM, 0.0d0, dbl_mb(k_eqme), 1)
c  Call EqmE to get field contributions due to QM e-
      call dimqm_EqmE2(rtdb, g_dens, geom, basis, 
     $         dbl_mb(k_eqme), dbl_mb(k_dimxyz), nDIM, g_eqme)
c      if(id.eq.0) write(luout,*)'back from eqme2'
c  Add eqme and eqmn and place into etot
      call ga_add(1.0d0, g_eqme, 1.0d0, g_eqmn, g_etot)
      call dimqm_getlfirst(rtdb, lfirst)
      if(.not.lfirst) then
         if(.not.rtdb_get(rtdb, 'dimqm:muind',mt_dbl,
     $                                    3*nDIM, dbl_mb(k_muold)))
     $      call errquit('main2: old dipoles not read',1,RTDB_ERR)
      else
         call dfill (3*nDIM, 0.0d0, dbl_mb(k_muold),1)
      endif
c
c      call dimqm_f2d(rtdb, Etot, muind, muold, dimxyz, 1, 'g', nDIM)
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
      call dimqm_iterative2(rtdb, g_etot, dbl_mb(k_dimxyz), nDIM)
c  ***** on exit, g_etot now holds muind *****
ccc Calculate interaction energy and put it in rtdb
      edimqm = -0.5d0 * ga_ddot(g_etot, g_eqmn)
      edimqm = edimqm + (-0.5d0 * ga_ddot(g_etot, g_eqme))
      if(.not.rtdb_put(rtdb, 'dimqm:edimqm',mt_dbl,1,edimqm))
     $   call errquit('main2: edimqm put failed',1,RTDB_ERR)
c  TODO: store ind. dip. in rtdb for now, so they can be accessed
c        when creating potential
c        ** NOTE ** it should be looked into for having the potential 
c        created and stored to rtdb instead. should be much smaller
c        than the dipoles since it only needs to apply to QM grid   
      call ga_init_fence()
      call ga_get(g_etot, 1, nOrder, 1, 1, dbl_mb(k_muind), 1)
      call ga_fence()
c      call ga_print(g_etot)
c      if (id.eq.0) then
c         do i = 0, nOrder-1
c            write(luout,*)'k_muind',i+1,' ',dbl_mb(k_muind+i)
c         enddo
c      endif
      if(.not.rtdb_put(rtdb, 'dimqm:muind',mt_dbl,nOrder,
     $   dbl_mb(k_muind))) call errquit('dim_main2: rtdb put muind
     $   failed',1,RTDB_ERR)
c  TODO: Check that dipoles are converged

c Place this cycle's dipoles into g_muold
      call ga_copy(g_etot, g_muold)

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
c      call dimqm_setlfirst(rtdb,.false.)
      stat = .true.
c      status = status.and.ma_free_heap(l_eqmn)
      stat = stat.and.ma_free_heap(l_eqme)
c      status = status.and.ma_free_heap(l_etot)
      stat = stat.and.ma_free_heap(l_muind)
c      status = status.and.ma_free_heap(l_muold)
      stat = stat.and.ma_chop_stack(l_dimxyz)
      stat = stat.and.ga_destroy(g_eqme)
      stat = stat.and.ga_destroy(g_eqmn)
      stat = stat.and.ga_destroy(g_etot)
      if (.not.stat) 
     $   call errquit('dimqm_main2:free stack failed',1,MA_ERR)

      if(lastcy) then
         if(id.eq.0) then
c  TODO
            write(luout,*)
     $'---------------------------- DIM/QM RESULTS -------------------'
            write(luout,*)''
            write(luout,*)''
            write(luout,*)
     $'-------------------------- DIM/QM Energy -----------------------'
            write(luout,*)'Total:   ', edimqm
            write(luout,*)'it seems the job finished, do this later'
         endif
c         stat = stat.and.ga_destroy(g_muind)
c         stat = stat.and.ga_destroy(g_eqmn)
c         stat = stat.and.ga_destroy(g_eqme)
c         stat = stat.and.ga_destroy(g_etot)
c         if (.not.stat) call errquit('main2: ga_dest failed',1,GA_ERR)
      endif
      end subroutine dimqm_main2


      subroutine dimqm_dip_err(muold, muind, nDIM, maxerr)
c-----------------------------------------------------------------------
c  Author: Jeff Becca
c  Calls: Nothing
c  Called from: dimqm_main2
c  Purpose: get max err from allocated arrays
      implicit none
c------------------------------ Input variables ------------------------
      double precision muold(3,nDIM)
      double precision muind(3,nDIM)
      integer nDIM
      double precision maxerr
      maxerr = MAXVAL(ABS(muold-muind))
      end subroutine dimqm_dip_err

      function dimqm_energy2(muind, eqme, eqmn, nDIM)
c-----------------------------------------------------------------------
c  Author: Jeff Becca
c  Calls: Nothing
c  Called from: dimqm_main2
c  Purpose: Calculate DIM/QM interaction energy
c-----------------------------------------------------------------------
      implicit none
c-------------------------- Input variables ----------------------------
      double precision muind(3*nDIM)
      double precision eqme(3*nDIM)
      double precision eqmn(3*nDIM)
      double precision dimqm_energy2, tmp1, tmp2
      integer nDIM, i

c      tmp1 = -0.5d0 * SUM(muind * eqme)
c      tmp2 = -0.5d0 * SUM(muind * eqmn)
c      dimqm_energy2 = tmp1 + tmp2
      tmp1 = 0.0d0
      tmp2 = 0.0d0
      do i = 1, 3*nDIM
         tmp1 = tmp1 + (muind(i)*eqme(i))
         tmp2 = tmp2 + (muind(i)*eqmn(i))
      enddo
      dimqm_energy2 = -0.5d0*(tmp1 + tmp2)
      return
      end function dimqm_energy2
