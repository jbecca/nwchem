      subroutine dimqm_energy(rtdb, nAtoms, muind, eqme, eqmn)
      implicit none
#include "stdio.fh"
#include "util.fh"
#include "dimqm_constants.fh"
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
c  jbecca: START
c  jbecca: END
      integer nAtoms
      double precision muind(3,nAtoms)
      double precision eqme(3, nAtoms)
      double precision eqmn(3, nAtoms)
c    Local Variables
      double precision en_dipole_el
      double precision en_dipole_nuc
      double precision efld_int
      double precision energy
      character*50 d
      character*50 dd
      integer rtdb
      logical lefield
      double precision finiteFieldX, finiteFieldY, finiteFieldZ
      
      if(ga_nodeid().ne.0) return
c
      dd =
     $ '=================================================='
      d =
     $ '--------------------------------------------------'
c
c    Calculate the DIM energy
      en_dipole_el = -HALF * SUM(muind * eqme)
      en_dipole_nuc= -HALF * SUM(muind * eqmn)
      efld_int  =  ZERO
      call dimqm_getlefield(rtdb, lefield)
      if (lefield) then
         if(.not.rtdb_get(rtdb,'dimqm:finiteFieldX',mt_dbl,
     $    1,finiteFieldX)) call errquit('get fieldX failed',1,RTDB_ERR)
         if(.not.rtdb_get(rtdb,'dimqm:finiteFieldY',mt_dbl,
     $    1,finiteFieldY)) call errquit('get fieldY failed',1,RTDB_ERR)
         if(.not.rtdb_get(rtdb,'dimqm:finiteFieldZ',mt_dbl,
     $    1,finiteFieldZ)) call errquit('get fieldZ failed',1,RTDB_ERR)
         efld_int  =  SUM(muind(1,:) * finiteFieldX)
     $             +  SUM(muind(2,:) * finiteFieldY)
     $             +  SUM(muind(3,:) * finiteFieldZ)
         efld_int  =  -HALF * efld_int
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

      end subroutine dimqm_energy
