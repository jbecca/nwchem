      subroutine dimqm_seed_init_tolerance(rtdb, request)
      implicit none
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
      integer rtdb
      double precision request
      logical dimqm_noseed
      logical ldebug
      logical dimqm_seeded
      double precision dimtol, dimtol0
      call dimqm_get_debug(rtdb,ldebug)
c   Initialize finding the seed
      dimqm_seeded = .false.
      call dimqm_put_seeded(rtdb, dimqm_seeded)
c   Return if seeding is off
      call dimqm_get_noseed(rtdb, dimqm_noseed)
      if(dimqm_noseed) return
c   If the user set tolerance is less than the pre-seed tolerance, save the user's
c   setting and change to the pre-seed tolerance.  This will be reversed once seeding
c   has started.
      if (.not. rtdb_get(rtdb, 'dimqm:dimtol', mt_dbl, 1, dimtol))
     $   call errquit('dim_seed: get dimtol failed',1,RTDB_ERR)
      if(dimtol < request) then
        dimtol0 = dimtol
        dimtol = request
        if (.not. rtdb_put(rtdb, 'dimqm:dimtol',mt_dbl, 1, dimtol))
     $   call errquit('dimqm_seed:put dimtol failed',1,RTDB_ERR)
        if (.not. rtdb_put(rtdb, 'dimqm:dimtol0',mt_dbl, 1, dimtol0))
     $   call errquit('dimqm_seed:put dimtol0 failed',1,RTDB_ERR)
        if(ga_nodeid().eq.0 .and. ldebug) then
          write(luout,*) "User requested tolerance below", request
          write(luout,*) "Setting pre-seed tolerance to", request
        end if
      end if
      end subroutine dimqm_seed_init_tolerance

      subroutine dimqm_check_dipoles(rtdb, tolerance, error)
      implicit none
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
      logical dimqm_noseed
      integer rtdb
      logical lfirst
      logical ldebug
      double precision tolerance, error, dimtol0, dimtol
      logical dimqm_seeded
      call dimqm_get_debug(rtdb,ldebug)
      call dimqm_getlfirst(rtdb,lfirst)
c      call dimqm_get_seeded(rtdb, dimqm_seeded)
c   Return if seeding is off
      if(dimqm_noseed) return
c   Return if seed already selected
      if(dimqm_seeded) return
c   Return if this is the first cycle
      if(lfirst) return
c   Check to see if the error is below our seeding tolerance,
      if(error < tolerance) then
c       Return if we've already set the seed
c       Otherwise, set the seed and revert the tolerance back to user requested level
        dimqm_seeded = .true.
   
        call dimqm_put_seeded(rtdb, dimqm_seeded)

        if(ga_nodeid().eq.0 .and. ldebug) then
          write(luout,*) "Dipole error below ", tolerance
          write(luout,*) "Setting current system as seed"
          write(luout,*) "Reverting tolerance back to", dimtol0
        end if
        dimtol = dimtol0
        if (.not. rtdb_put(rtdb, 'dimqm:dimtol', mt_dbl, 1, dimtol))
     $      call errquit('dim_check_dip: put dimtol failed',1,RTDB_ERR)
      end if
      end subroutine dimqm_check_dipoles

      subroutine dimqm_check_scf(rtdb, de, rms, derr, ipol, tolerance)
      implicit none
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"

      double precision de, rms(2), tolerance, derr(2), dimtol
      integer ipol
      logical e, d1, d2, g1, g2
      logical dimqm_noseed
      logical ldimqm
      integer rtdb
      logical ldebug
      logical dimqm_seeded
c   Return if dimqm is off

      call dimqm_get_debug(rtdb, ldebug)
      call dimqm_used(rtdb, ldimqm)
      call dimqm_get_seeded(rtdb, dimqm_seeded)
      if (.not.ldimqm) return
c   Return if seeding is off
      if(dimqm_noseed) return
c   Return if the seed hasn't been set yet
      if(.not.dimqm_seeded) return
      e = .false.
      d1 = .false.
      d2 = .false.
      g1 = .false.
      g2 = .false.
c   Check to see if the SCF RMS-density and delta E are below a given tolerance
c   If not, then this is not a good seed system and we reset the seed flag
      if(de .gt. tolerance) e = .true.
      if(dsqrt(rms(1)) .gt. tolerance) d1 = .true.
      if(derr(1) .gt. tolerance) g1 = .true.
      if(ipol .gt. 1) then
        if(dsqrt(rms(2)) .gt. tolerance) d2 = .true.
        if(derr(2) .gt. tolerance) g2 = .true.
      end if
      if(e .or. d1 .or. d2 .or. g1 .or. g2) then
        if(ga_nodeid().eq.0 .and. ldebug) then
          write(luout,*) "Seed was set, but some SCF error"
          write(luout,*) " is above ", tolerance
          write(luout,*) "Unsetting seed and reverting tolerance back"
          write(luout,*) "to 1.0d-4"
        end if
        dimqm_seeded = .false.
        call dimqm_put_seeded(rtdb, dimqm_seeded)
        dimtol = 1.0d-4
        if (.not. rtdb_put(rtdb, 'dimqm:dimtol', mt_dbl, 1, dimtol))
     $      call errquit('dim_check_dip: put dimtol failed',1,RTDB_ERR)
      end if
      end subroutine dimqm_check_scf
