      subroutine dimqm_input2(rtdb)
*
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "dimqm_constants.fh"
c
c     Input module for DIM/QM control parameters
c     ------------------------------------------
c          
c     Called from: input_parse.F
c
c     Calls: None
c
c     This subroutine reads the user defined parameters from the
c     dimqm block of the input file.  Defaults are set for
c     optional commands.  All global variables are defined and
c     stored in common blocks located within dimqm.fh
c
c     ----- Input -----
      integer rtdb ! Run time database
c
c     ----- Local Variables ----
c
      logical dimqm_noseed
      logical ldimqm
      logical lfirst
      integer nfield
      character*255 field
      character*4 scrnTmp
      logical lclfld
      logical lefield
      logical lrsp
      logical ldebug
      logical latmdip
      logical dimqm_go
      integer scrnType
      double precision dimtol
      double precision scrnFactor
      double precision finiteFieldX, finiteFieldY, finiteFieldZ
      integer id
c
c     ----- Defaults -----
c
      ldimqm            = .true.
      dimqm_go          = .false.
      ldebug            = .false. ! Print off DIM/QM debug info
      scrnType          = NOSCRN  ! DIM/QM screening function to use
      dimtol            = 1.0d-6  ! Iterative solver cutoff
      scrnFactor        = 1.0d0   ! DIM/QM screening factor
      lfirst            = .true.  ! This is the first SCF cycle
      dimqm_noseed      = .false. ! Seeding is on by default
      lrsp              = .true.
      lefield           = .false.
      latmdip           = .true.
c
c     ----- Read values from input -----
c
      id = ga_nodeid()
      if(id.eq.0) write(luout,*)'starting dimqm_input2'
 10   if(.not.inp_read())
     $     call errquit('dimqm_input: unexpected eof',911, INPUT_ERR)
      nfield = inp_n_field()
 20   if (.not. inp_a(field))
     $     call errquit('dimqm_input: failed to read field',911,
     &       INPUT_ERR)
c
      if(inp_compare(.false.,'end',field)) then
         go to 40
c
c   --- Methods
c
      else if(inp_compare(.false.,'off', field)) then
         ldimqm = .false.
      else if(inp_compare(.false.,'tolerance',field)) then
         if(.not. inp_f(dimtol))
     $      call errquit('dimqm_input: tol',911, INPUT_ERR)
c
c   --- Interaction Options
c
      else if(inp_compare(.false.,'efield',field)) then
         lefield = .true.
         if(.not. inp_f(finiteFieldX))
     $     call errquit('dimqm_input: efieldX val', 1, INPUT_ERR)
         if(.not. inp_f(finiteFieldY))
     $     call errquit('dimqm_input: efieldY val', 1, INPUT_ERR)
         if(.not. inp_f(finiteFieldZ))
     $     call errquit('dimqm_input: efieldZ val', 1, INPUT_ERR)

      else if(inp_compare(.false.,'screen',field)) then
         if(.not. inp_a_trunc(scrnTmp))
     $      call errquit('dimqm_input: screen type',911, INPUT_ERR)
         if(scrnTmp .eq. "NONE")  scrnType = NOSCRN
         if(scrnTmp .eq. "EXP ")  then
             scrnType = EXPSCRN
             if(.not. inp_f(scrnFactor))
     $         call errquit('dimqm_input: screen factor',911, INPUT_ERR)
         end if
         if( scrnTmp .eq. "ERF " )  then
           scrnType = ERFSCRN
           if( .not. inp_f(scrnFactor) )
     $       call errquit('dimqm_input: screen factor',911, INPUT_ERR)
         end if
c      Seeding
      else if(inp_compare(.false.,'noseed',field)) then
         dimqm_noseed = .true.
c
      else if(inp_compare(.false.,'localfield',field)) then
         lclfld = .true.
      else if(inp_compare(.false.,'noresp',field)) then
         lrsp = .false.
c   --- Printing Options
c
      else if(inp_compare(.false.,'debug',field)) then
         ldebug = .true.
      else if(inp_compare(.false.,'printatomicdipoles',field)) then
         latmdip = .true.
      endif
      go to 10
c
 40   continue

c     Some DIM/QM global variables are needed in the rtdb to
c     reduce module dependency contamination
      if (.not. rtdb_put(rtdb, 'dimqm:ldimqm', mt_log, 1, ldimqm))
     $  call errquit('dimqm_input: ldimqm put failed', 1, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'dimqm:dimqm_go', mt_log, 1, dimqm_go))
     $  call errquit('dimqm_input: dimqm_go put failed', 1, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'dimqm:ldebug', mt_log, 1, ldebug))
     $  call errquit('dimqm_input: ldebug put failed', 1, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'dimqm:scrnType', mt_int, 1, scrnType))
     $  call errquit('dimqm_input: scrnType put failed', 1, RTDB_ERR)
      if (.not. rtdb_put(rtdb,'dimqm:scrnFactor',mt_dbl,1,scrnFactor))
     $  call errquit('dimqm_input: scrnFactor put failed', 1, RTDB_ERR)
c -------------------------------------------------------------------
      if (.not. rtdb_put(rtdb, 'dimqm:lfirst', mt_log, 1, lfirst))
     $  call errquit('dimqm_input: lfirst put failed', 1, RTDB_ERR)
      if(.not.rtdb_put(rtdb,'dimqm:noseed',mt_log,1,dimqm_noseed))
     $  call errquit('dimqm_input: noseed put failed', 1, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'dimqm:dimtol', mt_dbl, 1, dimtol))
     $  call errquit('dimqm_input: dimtol put failed', 1, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'dimqm:lclfld', mt_log, 1, lclfld))
     $  call errquit('dimqm_input: lclfld put failed', 1, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'dimqm:lrsp', mt_log, 1, lrsp))
     $  call errquit('dimqm_input: lrsp put failed', 1, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'dimqm:lefield', mt_log, 1, lefield))
     $  call errquit('dimqm_input: lefield put failed', 1, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'dimqm:latmdip', mt_log, 1, latmdip))
     $  call errquit('dimqm_input: latmdip put failed', 1, RTDB_ERR)
      return
      end subroutine dimqm_input2
