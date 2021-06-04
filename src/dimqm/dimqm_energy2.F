      subroutine dimqm_energy2(rtdb, g_muind, g_eqme, g_eqmn)
      implicit none
#include "stdio.fh"
#include "util.fh"
#include "dimqm_constants.fh"
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "dimqm.fh"

c     Input variables
      integer rtdb
      integer g_muind
      integer g_eqme
      integer g_eqmn

c    Local Variables
      double precision en_dipole_el
      double precision en_dipole_nuc
      double precision efld_int
      double precision energy
      character*50 d
      character*50 dd
      logical lefield, stat
      double precision finiteFieldX, finiteFieldY, finiteFieldZ
      integer l_fld, k_fld
      integer g_efield
      
      if(ga_nodeid().ne.0) return
c
      dd =
     $ '=================================================='
      d =
     $ '--------------------------------------------------'
c
c    Calculate the DIM energy
      efld_int  =  ZERO
      en_dipole_nuc = -0.5d0 * ga_ddot(g_DIMmuind, g_eqmn)
      en_dipole_el  = -0.5d0 * ga_ddot(g_DIMmuind, g_eqme)

      call dimqm_getlefield(rtdb, lefield)
      if (lefield) then
         if(.not.ma_alloc_get(mt_dbl,3*DIM_nAtoms,'dimqm:fld',l_fld,
     $   k_fld)) call errquit('dimqm_gs: malloc k_fld failed',
     $   1,MA_ERR)
        if(.not.ga_create(mt_dbl,DIM_nOrder,1,'g_efield',0,0,g_efield))
     $   call errquit('dimqm_gs: create g_efield failed',0,GA_ERR)
         if(.not.rtdb_get(rtdb,'dimqm:finiteFieldX',mt_dbl,
     $    1,finiteFieldX)) call errquit('get fieldX failed',1,RTDB_ERR)
         if(.not.rtdb_get(rtdb,'dimqm:finiteFieldY',mt_dbl,
     $    1,finiteFieldY)) call errquit('get fieldY failed',1,RTDB_ERR)
         if(.not.rtdb_get(rtdb,'dimqm:finiteFieldZ',mt_dbl,
     $    1,finiteFieldZ)) call errquit('get fieldZ failed',1,RTDB_ERR)

         call dfill(DIM_nAtoms, finiteFieldX, dbl_mb(k_fld), 3) 
         call dfill(DIM_nAtoms, finiteFieldY, dbl_mb(k_fld+1), 3) 
         call dfill(DIM_nAtoms, finiteFieldZ, dbl_mb(k_fld+2), 3) 
         call ga_put(g_efield,1,3*DIM_nAtoms,1,1,dbl_mb(k_fld), 1)

         efld_int  =  -0.5d0 * ga_ddot(g_DIMmuind, g_efield)
         
         stat = .true.
         stat = stat.and.ga_destroy(g_efield)
         stat = stat.and.ma_free_heap(l_fld)
         if(.not.stat) 
     $      call errquit('dimqm_energy: efield clean failed',1,GA_ERR)

      endif
      energy = en_dipole_el + en_dipole_nuc
c
c
      write(luout,*) dd
      write(LuOut,*) 'DIM/QM Energy'
      write(LuOut,*) d
      write(LuOut,111) en_dipole_el
      write(LuOut,222) en_dipole_nuc
c      write(LuOut,444) efld_int
      write(LuOut,333) energy

      write(luout,*) dd
      write(luout,*) ''
      if (lefield) then 
         write(luout,*) dd
         write(luout,*)'DIM System Energy'
         write(LuOut,*) d
         write(luout,555) efld_int
         write(luout,*) dd
      endif

      call util_flush(LuOut)
 111  format(' Polarization (Dipole, el.) = ',f20.10)
 222  format(' Polarization (Dipole, nuc) = ',f20.10)
 333  format('              Total Energy  = ',f20.10)
 444  format(' Interaction with Efield    = ',f20.10)
 555  format('                      Total = ',f20.10)

      end subroutine dimqm_energy2
