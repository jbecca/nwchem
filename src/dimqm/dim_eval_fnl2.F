      subroutine dim_eval_fnl2(ptnl, nq, qwght, xyz, qxyz, 
     $                           muind, nDIM, rtdb)
      implicit none
c-----------------------------------------------------------------------
c  Purpose: Calculate final g.s. DIM potential at each quad. point
c  Calls: none
c  Called by: dim_grid_quadv0b
c-----------------------------------------------------------------------
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"

c  Input Variables
      double precision ptnl(nq)
      double precision qwght(nq)
      double precision xyz(3, nDIM)
      double precision qxyz(3, nq)
      double precision muind(3, nDIM)
      integer nq, nDIM, rtdb

c  Local Variables
      integer n, m, id
      double precision r(3), dist, dist3

      id = ga_nodeid()

      call dfill(nq, 0.0d0, ptnl, 1)

      do n = 1, nq
         do m = 1, nDIM
         r(:) = xyz(:,m) - qxyz(:,n) 
         dist = SQRT((r(1)**2 + r(2)**2 + r(3)**2))
         if (dist < 1.0d-12) dist = 1.0d-12
         dist3 = dist**3
         ptnl(n) = ptnl(n) - muind(1,m)*r(1)/dist3
         ptnl(n) = ptnl(n) - muind(2,m)*r(2)/dist3
         ptnl(n) = ptnl(n) - muind(3,m)*r(3)/dist3
         enddo          !m
      enddo             !n

c  combine with grid weights
      do n = 1, nq
         ptnl(n) = -ptnl(n)*qwght(n)
      enddo             !n 

c      if (id.eq.0) write(luout,*)'ptnl',ptnl
      return
      end subroutine dim_eval_fnl2
