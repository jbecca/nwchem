      subroutine dimqm_EqmN2(rtdb, geom, xyz, g_eqmn, nDIM)
c-----------------------------------------------------------------------
c  Author: Jeff Becca
c  Calls: dimqm_fldpot
c  Called from: dimqm_main2
c  Purpose: get field contribution from QM nuclei
c-----------------------------------------------------------------------
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "msgids.fh"
c------------------------------Input variables--------------------------
      integer rtdb, geom, nDIM, g_eqmn
      double precision xyz(3,nDIM)
c      double precision eqmn(3,nDIM)
c-------------------------------Local variables-------------------------
      integer id, nNuc, l_qmxyz, k_qmxyz, l_qmchrg, k_qmchrg
      integer mLo, mUp, mOnNode, m, n
      logical ldebug

      id = ga_nodeid()
      ldebug = .true.

c  gather QM geom data
      if(.not.geom_ncent(geom, nNuc))
     $   call errquit('dimqm_eqmn2: geom ncent?',1,GEOM_ERR)

c  allocate memory for qm coords
      if(.not.ma_push_get(mt_dbl,3*nNuc,'QM xyz',l_qmxyz,k_qmxyz))
     $   call errquit('EqmN2: malloc k_qmxyz failed',1,MA_ERR)
      if(.not.ma_push_get(mt_dbl,nNuc,'QM chrg',l_qmchrg,k_qmchrg))
     $   call errquit('EqmN2: malloc k_qmchrg failed',1,MA_ERR)

c  get values for arrays
      if(.not.geom_cart_coords_get(geom, dbl_mb(k_qmxyz)))
     $   call errquit('EqmN2: get coords failed',1,GEOM_ERR)
      if(.not.geom_cart_get_charges(geom, nNuc, dbl_mb(k_qmchrg)))
     $   call errquit('EqmN2: get charges failed',1,GEOM_ERR)

      call dimqm_fldpot2(rtdb, xyz, dbl_mb(k_qmxyz), 
     $                        dbl_mb(k_qmchrg),g_eqmn, nNuc, nDIM)
c  TODO find a way to use all GA's and use ga operations
c  Deallocate mem
      if(.not.ma_chop_stack(l_qmxyz))
     $   call errquit('EqmN2: chop l_qmxyz failed',1,MA_ERR)

      end subroutine dimqm_EqmN2

      

