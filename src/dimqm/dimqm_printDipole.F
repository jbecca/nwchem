      subroutine dimqm_printDipole(rtdb, muind, lpprint, nDIM)
c
c
c
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

c   Input variables
      integer rtdb
      integer g_muind
c      double precision qind(nDIM, icmplx)
      logical lpprint
      integer nDIM
      double precision muind(3, nDIM)
c
c   Local variables
      double precision diptot(3)
      character*60 dd
      character*60 d

      if(ga_nodeid().ne.0) return

      diptot = SUM(muind, DIM = 2)
      dd =
     $ '============================================================='
      d =
     $ '-------------------------------------------------------------'
c   BLAS matrix-matrix multiply to add in charge term

      if(lpprint) then
       write(luout,'(/1x,a)') dd
       write(luout,'(1x,a)')
     $                   'Total induced dipole moment in DIM system :'
       write(luout,'(1x,a)') d
       write(luout,'(2x,a,2x,3f16.8)') 'A.U.:', diptot(:)
       write(luout,'(2x,a,1x,3f16.8)') 'Debye:', diptot(:)*AU2DEBYE
       write(luout,'(1x,a/)') dd
      else
       write(luout,'(1x,a)') 
     $                   'Total induced dipole moment in DIM system :'
       write(luout,'(3f14.8)') diptot(:)
      end if
      
      end subroutine dimqm_printDipole
