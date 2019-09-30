      subroutine dimqm_direct(rtdb, dimxyz, g_fld, nDIM, nDIMTypes)
c----------------------------------------------------------------------
c  Author: Jeff Becca
c  Calls: 
c  Called from: dimqm_f2d2 (in the future)
c  Purpose: Constructs the full interaction matrix for the DIM system
c           and calls routine for Cholesky decomp.
c-----------------------------------------------------------------------
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
c------------------------- Input Variables -----------------------------
      integer rtdb, g_fld, g_muind, nDIM, nDIMTypes
c------------------------- Local Variables -----------------------------
      integer id, tx(nDIM), nOrder, n, m, m1, m2, m3, n1, n2, n3
      logical ldebug
      double precision sPol(nDIMTypes)
      double precision bound(nDIMTypes)
      double precision dimxyz(3, nDIM)
      double precision drudeterms(2, nDIMTypes)
      double precision inv_p, dist, invdist
      double precision r(3)
      logical hasExpDie(nDIMTypes)
      double precision t2(3,3)
      integer g_int_mat, final_prod
      interface
        function t2r(r, dist, invdist, a1, a2, pol1, pol2)
          integer, intent(in) :: a1, a2
          double precision, intent(in) :: dist, invdist
          double precision, intent(in) :: pol1, pol2
          double precision t2r(3,3)
          double precision, intent(in) :: r(3)
        end function t2r
      end interface

      id = ga_nodeid()
      nOrder = 3*nDIM

      if(.not. rtdb_get(rtdb,'dimqm:spol', mt_dbl, nDIMTypes, sPol))
     $  call errquit('dimqm:spol rtdb get failed', 0, RTDB_ERR)
      if(.not. rtdb_get(rtdb,'dimqm:typeindex', mt_int, nDIM, tx))
     $  call errquit('dimpar:type index rtdb get failed', 0, RTDB_ERR)
      if(.not. rtdb_get(rtdb,'dimqm:bound', mt_dbl, nDIMTypes, bound))
     $  call errquit('dimqm:bound rtdb get failed', 0, RTDB_ERR)
      if(.not. rtdb_get(rtdb,'dimqm:drude', mt_dbl, nDIMTypes*2,
     $                                    drudeterms))
     $  call errquit('dimqm:drude rtdb get failed', 0, RTDB_ERR)
      if(.not. rtdb_get(rtdb,'dimqm:hasexp', mt_log, nDIMTypes,
     $                                    hasExpDie))
     $  call errquit('dimqm:hasexp rtdb get failed', 0, RTDB_ERR)

      if(.not.ga_create(mt_dbl, nOrder, nOrder, 'dimqm:int_mat', 
     $   nOrder, -1, g_int_mat)) call errquit(
     $   'dimqm_direct: create g_int_mat failed',0,GA_ERR)
      call ga_zero(g_int_mat)

c  This all only works for static case as of now, no attention
c  is paid to the frequency dependence
      do m = 1, nDIM
         m1 = 3*m-2
         m2 = 3*m-1
         m3 = 3*m
c  Get the diagonal elements first
         inv_p = 1.0d0 / sPol(tx(m))
         call ga_put(g_int_mat, m1, m1, m1, m1, inv_p, 1)
         call ga_put(g_int_mat, m2, m2, m2, m2, inv_p, 1)
         call ga_put(g_int_mat, m3, m3, m3, m3, inv_p, 1)
c  now get off diagonal elements
         do n = 1, nDIM
            n1 = 3*n-2
            n2 = 3*n-1
            n3 = 3*n
c  cycle if diagonal
            if (m == n) cycle 
            r(1:3) = dimxyz(1:3,n) - dimxyz(1:3,m)
            dist = SQRT(r(1)**2 + r(2)**2 + r(3)**2)
            invdist = 1.0d0 / dist
            t2 = t2r(r, dist, invdist, m, n, sPol(tx(m)), sPol(tx(n)))
            call ga_put(g_int_mat, m1, m3, n1, n3, t2, 1) 
         enddo       !n
      enddo          !m
      call ga_sync
      final_prod = ga_solve(g_int_mat, g_fld)
      if (final_prod.ne.0) then
         call errquit('dim_direct: solver failed',0,GA_ERR)
      endif
      call ga_sync
      end subroutine dimqm_direct
