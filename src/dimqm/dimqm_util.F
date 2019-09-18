c     
c     $Id: dft_scf.F 28382 2016-07-21 01:23:33Z niri $
c     
      subroutine dimqm_used(rtdb, ldimqm)
      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      logical ldimqm
      integer rtdb
      if (.not.rtdb_get(rtdb, 'dimqm:ldimqm',mt_log,1,ldimqm)) then
         ldimqm = .false.
      endif
      return
      end
      subroutine dimqm_getgo(rtdb, dimqm_go)
      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      logical dimqm_go
      integer rtdb
      if (.not.(rtdb_get(rtdb, 'dimqm:dimqm_go',mt_log,1,dimqm_go)))
     $      call errquit('dimqm_getgo failed',1,RTDB_ERR)
      return
      end
      subroutine dimqm_setgo(rtdb, dimqm_go)
      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      logical dimqm_go
      if (.not.(rtdb_put(rtdb, 'dimqm:dimqm_go',mt_log,1,dimqm_go)))
     $      call errquit('dimqm_setgo failed',1,RTDB_ERR)
      return
      end
      subroutine dimqm_setlfirst(rtdb,lfirst)
      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      logical lfirst
      if (.not.(rtdb_put(rtdb, 'dimqm:lfirst',mt_log,1,lfirst)))
     $      call errquit('dimqm_setlfirst failed',1,RTDB_ERR)
      return
      end
      subroutine dimqm_getlrsp(rtdb,lrsp)
      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      logical lrsp
      if (.not.(rtdb_get(rtdb, 'dimqm:lrsp',mt_log,1,lrsp)))
     $      call errquit('dimqm_getlrsp failed',1,RTDB_ERR)
      return
      end
      subroutine dimqm_getlcpks(rtdb,lcpks)
      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      logical lcpks
      if (.not.(rtdb_get(rtdb, 'dimqm:lcpks',mt_log,1,lcpks)))
     $      lcpks = .false.
      return
      end
      subroutine dimqm_getlclfld(rtdb,lclfld)
      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      logical lclfld
      if (.not.(rtdb_get(rtdb, 'dimqm:lclfld',mt_log,1,lclfld)))
     $      call errquit('dimqm_getlclfld failed',1,RTDB_ERR)
      return
      end subroutine dimqm_getlclfld

c      subroutine dimqm_init(rtdb)
c      implicit none
c#include "dimqm.fh"
c#include "rtdb.fh"
c#include "mafdecls.fh"
c      integer rtdb ! Run time database
cc     Some DIM/QM global variables are needed in the rtdb to
c      if (.not. rtdb_get(rtdb, 'dimqm:ldimqm', mt_log, 1, ldimqm))
c     A     ldimqm=.false.
c      if (.not. rtdb_get(rtdb, 'dimqm:lrsp', mt_log, 1, lrsp))
c     A     lrsp=.false.
c      if (.not. rtdb_get(rtdb, 'dimqm:lfirst', mt_log, 1, lfirst))
c     A     lfirst=.false.
c      if (.not. rtdb_get(rtdb, 'dimpar:ndim', mt_int, 1, ndim))
c     A     ndim=0
c      if (.not. rtdb_get(rtdb,'dimpar:ndimtypes', mt_int, 1, ndimtypes))
c     A     ndimtypes=0
c      if (.not. rtdb_get(rtdb, 'dimpar:ndim', mt_int, 1, ndim))
c     A     ndim=0
c      if (.not. rtdb_get(rtdb,'dimpar:norder', mt_int, 1, norder))
c     A     norder=0
c      if (.not. rtdb_get(rtdb,'dimqm:alg', mt_int, 1, alg))
c     A     alg=0
c      if (.not. rtdb_get(rtdb, 'dimqm:dimtol', mt_dbl, 1, dimtol))
c     $     dimtol=0d0
c      write(6,*)'ldimqm from dimqm_init',ldimqm
c      return
c      end

c The following subroutines are additions to things that need to be
c  set or read during DIM/QM calculations. They will be made individually
c  and also one master routine that can pull and put all relevant values

      subroutine dimqm_get_noseed(rtdb, dimqm_noseed)
      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      logical dimqm_noseed
      
      if (.not.(rtdb_get(rtdb, 'dimqm:noseed',mt_log,1,dimqm_noseed)))
     $      call errquit('dim_get_noseed failed',1,RTDB_ERR)
      return
      end subroutine dimqm_get_noseed

      subroutine dimqm_put_noseed(rtdb, dimqm_noseed)
      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      logical dimqm_noseed
      
      if (.not.(rtdb_put(rtdb,'dimqm:noseed',mt_log,1,dimqm_noseed)))
     $      call errquit('dim_get_noseed failed',1,RTDB_ERR)
      return
      end subroutine dimqm_put_noseed

      subroutine dimqm_getenergy(rtdb, edimqm)
      implicit none
      double precision edimqm
      integer rtdb
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      if (.not.(rtdb_get(rtdb, 'dimqm:edimqm',mt_dbl,1,edimqm)))
     $      call errquit('dim_getenergy failed',1,RTDB_ERR)
      return
      end subroutine dimqm_getenergy

      subroutine dimqm_putenergy(rtdb, edimqm)
      implicit none
      double precision edimqm
      integer rtdb
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      if (.not.(rtdb_put(rtdb, 'dimqm:edimqm',mt_dbl,1,edimqm)))
     $      call errquit('dim_putfailed',1,RTDB_ERR)
      return
      end subroutine dimqm_putenergy

      subroutine dimqm_get_debug(rtdb, ldebug)
      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      logical ldebug
      
      if (.not.(rtdb_get(rtdb, 'dimqm:ldebug',mt_log,1,ldebug)))
     $      call errquit('dim_get_debug failed',1,RTDB_ERR)
      return
      end subroutine dimqm_get_debug

      subroutine dimqm_getndim(rtdb, nDIM)
      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      integer nDIM
      
      if (.not.(rtdb_get(rtdb, 'dimpar:ndim',mt_int,1,nDIM)))
     $      call errquit('dimqm_getndim failed',1,RTDB_ERR)
      return
      end subroutine dimqm_getndim

      subroutine dimqm_getlfirst(rtdb,lfirst)
      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      logical lfirst
      if (.not.(rtdb_get(rtdb, 'dimqm:lfirst',mt_log,1,lfirst)))
     $      call errquit('dimqm_getlfirst failed',1,RTDB_ERR)
      return
      end subroutine dimqm_getlfirst

      subroutine dimqm_getlefield(rtdb,lefield)
      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      logical lefield
      if (.not.(rtdb_get(rtdb, 'dimqm:lefield',mt_log,1,lefield)))
     $      call errquit('dimqm_getlefield failed',1,RTDB_ERR)
      return
      end subroutine dimqm_getlefield

      subroutine dimqm_getlatmdip(rtdb,latmdip)
      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      logical latmdip
      if (.not.(rtdb_get(rtdb, 'dimqm:latmdip',mt_log,1,latmdip)))
     $      call errquit('dimqm_getlatmdip failed',1,RTDB_ERR)
      return
      end subroutine dimqm_getlatmdip

      subroutine dimqm_get_seeded(rtdb,dimqm_seeded)
      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      logical dimqm_seeded
      if (.not.(rtdb_get(rtdb, 'dimqm:seeded',mt_log,1,dimqm_seeded))) 
     $   call errquit('dimqm_getseeded failed',1,RTDB_ERR)
      return
      end subroutine dimqm_get_seeded

      subroutine dimqm_put_seeded(rtdb, dimqm_seeded)
      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      integer rtdb
      logical dimqm_seeded
      
      if (.not.(rtdb_put(rtdb,'dimqm:seeded',mt_log,1,dimqm_seeded)))
     $      call errquit('dimqm_put_seeded failed',1,RTDB_ERR)
      return
      end subroutine dimqm_put_seeded
