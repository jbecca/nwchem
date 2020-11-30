      subroutine dim_grid_gs(rtdb, g_vdim)
c
      implicit none
      integer rtdb
      integer g_vdim

      call dim_grid_gs_gen(rtdb, g_vdim, 1, 1, .false., 0, 1)

      return
      end
c
      subroutine dim_grid_gs_gen(rtdb, g_vdim,
     &      nmat, calc_type, l3d, ipm, imag)

      implicit none

#include "errquit.fh"
c
      integer rtdb
      integer g_vdim

      Logical oprint,oprintw
      Logical xc_chktau 
      external xc_rep_close, xc_chktau
      logical xc_rep_close
      integer nmat              ! [input] number of DIM matrices to make;
c                               !         only used for lhs cpks
      integer calc_type         ! what to calculate [input]
c                               !   1 = XC matrices for SCF
c                               !   2 = XC matrices for CPKS LHS
      logical l3d               ! true if 3d GA is used
      integer ipm               ! Control over dipoles
      integer imag              ! Control over dipoles
c
#include "cgridfile.fh"
#include "geom.fh"
#include "cdft.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "rtdb.fh"
#include "dft2drv.fh"
#include "steric.fh"
#include "dimqm.fh"
c
      logical rtdb_status       
      logical grid_reopen
      external grid_reopen
      logical lsa ! Qin
c
      integer me, nproc, nxyz, nctrs, 
     &        nbf_ao_mxnbf_ce2,  
     &        ncoef_max,
     &        npol
c
      integer lxyz, ixyz, lcharge, icharge, 
     &        lqxyz, iqxyz, lqwght, iqwght, 
     &        lniz, iniz, lexpo, iexpo, lcoef, icoef, lFmat, iFmat, 
     &        lPmat, iPmat, lrho, irho, ldelrho, idelrho, lBmat, iBmat,
     &        lhess, ihess, lrq, irq, 
     &        lfin, ifin, ltags, itags,
     &        ltau, itau,  lmmat, immat, lttau, ittau
c
      integer lcntoce, icntoce, lcntobfr, icntobfr,
     &        lcetobfr, icetobfr, lrdens_atom, irdens_atom,
     &        nscr, lscr, iscr,lamat,iamat,icmat,lcmat

      double precision nelec
c
      integer k_muind, l_muind
      integer libf_ao,iibf_ao,irqbuf,lrqbuf,lbp,ibp,
     .     irchi_atom,irdelchi_atom,lrchi_atom,lrdelchi_atom
      integer iqs,lqs,nqlen,mxnrad,
     , laqbuf, iaqbuf,lRij,iRij,ndbl,l_qat,k_qat,irqlen
      integer l_Amat2, k_Amat2, l_Cmat2, k_Cmat2, l_Mmat2, k_Mmat2
      integer len_rho, len_delrho, len_Amat2, len_Cmat2, len_Mmat2
      integer len_tau, len_mmat, len_ttau
      logical do_scf, do_cpks_l, do_cpks_r, do_nmrks_r
      integer max_elem_ang_scr, max_component
      integer nexpo,nibf
      double precision job_grid_acc,rho_err
      integer ii,isp
c
      logical fon ! fractional occupancy
      double precision pstrace  ! a measure of the total charge
      integer i3, idir
      integer l_dimxyz, k_dimxyz
      character*(1) direction(3)
      character*(1) dpm(2)
      character*(1) dri(2)
      data direction /'x', 'y', 'z'/
      data dpm /'+', '-'/
      data dri /'r', 'i'/
      integer nDIM
      logical ldebug

c
c      call dimqm_get_debug(rtdb, ldebug)
      ldebug = .true.
      me = ga_nodeid()
      nproc = ga_nnodes()
      oprint = ldebug 
c      if(oprint) then
c        write(luout,*) "dim_grid_quadv0 start"
c        write(*,*) "GEOM:", geom
c      end if
c      oprint = .true.
      oprintw = util_print('quadrature warning', print_medium)
c
      irdelchi_atom = 0
c
      npol = (ipol*(ipol+1))/2
      if (.not.rtdb_get(rtdb,'dft:largenode', mt_log, 1, largenode))
     &     largenode=.false.
c
c     Open grid pts file
c     
cng
      if (grid_written .and. n_rec_in_file .le. 0 ) then
c
c       Somebody (fock_xc?) closed the grid file, better reopen it
c       before continuing. Note, fock_xc tends to open the grid file
c       at the beginning and close it at the end.
c
        if (.not.grid_reopen(geom)) call errquit(
     &    'dim_grid_quadv0: could not open grid file',0,DISK_ERR)
      endif
cng
      if (.not. grid_written .or. n_per_rec .le. 0 ) then
         call grid_packing_info
      end if
c
c     Set up logicals denoting calculation type
c
      do_scf = calc_type .eq. 1
      do_cpks_l = calc_type .eq. 2
      if (.not. (do_scf .or. do_cpks_l)) then
         call errquit('dim_grid_quadv0_gen: illegal calculation type',
     &    1,INPUT_ERR)
      endif
c
c      if (oprint) then
c        write(*,*) " *** dim_grid_quadv0: do_scf:",do_scf 
c        write(*,*) " *** dim_grid_quadv0: n_per_rec:",n_per_rec
c      end if
c
      len_tau = 1
      len_mmat = 1
      len_ttau = 1
      len_Amat2 = 1
      len_Cmat2 = 1
      len_Mmat2 = 1
c
      if (do_scf) then
         len_rho = npol*n_per_rec
      elseif (do_cpks_l) then
         len_rho = n_per_rec*(nmat+1)+n_per_rec
      endif
      if (.not. geom_ncent(geom, nctrs)) 
     &   call errquit('grid_den: geom_ncent failed',geom, GEOM_ERR)
      nxyz = 3*nctrs
      if (.not.MA_Alloc_Get(MT_Dbl,nxyz,'xyz',lxyz,ixyz))
     &   call errquit('grid_den: cannot allocate xyz',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,nctrs,'charge',lcharge,icharge))
     &   call errquit('grid_den: cannot allocate charge',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Byte,nctrs*16,'tags',ltags,itags))
     &   call errquit('grid_den: cannot allocate tags',0, MA_ERR)
      if (.not. geom_cart_get(geom, nctrs, Byte_MB(itags),
     &                        Dbl_MB(ixyz), Dbl_MB(icharge)))
     &   call errquit('grid_den: geom_cart_get failed',74, GEOM_ERR)
c
c      write(*,*)' Coordinates from grid_den.F '
c      call output(dbl_mb(ixyz), 1, 3, 1, nctrs, 3, nctrs, 1)
c     &                    
c
      if (.not.ma_pop_stack(ltags))
     &   call errquit('grid_den: cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(lcharge))
     &   call errquit('grid_den: cannot pop stack',0, MA_ERR)
c
      nbf_ao_mxnbf_ce2 = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
c      if (do_nmrks_r) nbf_ao_mxnbf_ce2=nbf_ao_mxnbf_ce2*3
c
c     Divide the multi-center integrand into a set of single center integrands.
c
      if (.not.MA_Push_get(mt_dbl,3*n_per_rec,'qxyz',lqxyz,iqxyz))
     &   call errquit('grid_den: cannot allocate qxyz',0, MA_ERR)
      if (.not.MA_Push_get(mt_dbl,n_per_rec,'qwght',lqwght,iqwght))
     &   call errquit('grid_den: cannot allocate qwght',0, MA_ERR)
      if (.not.MA_Push_get(MT_int,nctrs,'atom list',lniz,iniz))
     &   call errquit('grid_den: cannot allocate atom list',0, MA_ERR)
      if (.not.MA_Push_get(MT_int,nctrs,'atom nbf',lfin,ifin))
     &   call errquit('grid_den: cannot allocate atom nbf',0, MA_ERR)
      nexpo=nbf_ao_mxprim
      if (.not.MA_Push_Get(MT_Dbl,nbf_ao_mxprim,'expo',lexpo,iexpo))
     &   call errquit('grid_den: cannot allocate expo',0, MA_ERR)
c
      ncoef_max = nbf_ao_mxnbf_cn * nbf_ao_mxprim
c
      if (.not.MA_Push_Get(MT_Dbl,ncoef_max,'coef',lcoef,icoef))
     &   call errquit('grid_den: cannot allocate coef',0, MA_ERR)
      if (.not.MA_Push_Get(MT_dbl,nbf_ao_mxnbf_ce2,'DENSF',lFmat,iFmat))
     &   call errquit('grid_den: cannot allocate DENSF',0, MA_ERR)
      if (.not.MA_Push_Get(MT_dbl,nbf_ao_mxnbf_ce2,'DENS',lPmat,iPmat))
     &   call errquit('grid_den: cannot allocate DENS',0, MA_ERR)
c
      if (.not.MA_Push_Get(MT_Dbl,len_rho,'rho',lrho,irho))
     &   call errquit('grid_den: cannot allocate rho',0, MA_ERR)
c
c
      if (.not.MA_Push_Get(mt_dbl,nbf_ao_mxnbf_ce*n_per_rec,'Bmat',
     &   lBmat,iBmat))call errquit('grid_den: cannot allocate Bmat',0,
     &       MA_ERR)
c
c     Assume this numerical density routine can be called from anywhere.
c     -compute mapping vectors
c 
      if (.not.ma_push_get
     &   (mt_int,nctrs*2,'cntoce map',lcetobfr,icetobfr))
     &   call errquit('grid_den:push_get failed', 13, MA_ERR)
      if (.not.ma_push_get
     &   (mt_int,nshells_ao,'cntoce map',lcntoce,icntoce))
     &   call errquit('grid_den:push_get failed', 13, MA_ERR)
      if (.not.ma_push_get
     &   (mt_int,nshells_ao*2,'cntoce map',lcntobfr,icntobfr))
     &   call errquit('grid_den:push_get failed', 13, MA_ERR)
c     
      call build_maps(ao_bas_han, int_mb(icntoce), int_mb(icntobfr), 
     &                int_mb(icetobfr), nctrs, nshells_ao)
      if (.not.ma_chop_stack(lcntoce))
     &   call errquit('grid_den: cannot pop stack',1, MA_ERR)
c
c     allocate scratch array for grid compression
c
      nibf=nbf_ao
      if (.not.MA_Push_Get(mt_int, nibf, 'ibf_ao', libf_ao,
     &     iibf_ao))
     &     call errquit('grid_quadv0: cannot allocate ibf_ao',2,
     &       MA_ERR)
      if (.not.MA_Push_get(MT_dbl, 4*buffer_size+4, 
     &                     'quad pts buffer', lrqbuf, irqbuf))
     &   call errquit('grid_gen: cannot allocate quad pt buffer', 3,
     &       MA_ERR)

      if (.not.MA_Push_Get(MT_Dbl,nctrs,'rchi_atom',
     &     lrchi_atom,irchi_atom))
     &     call errquit('dft_scf: cannot allocate rchi_atom',0, MA_ERR)
      irqlen=nctrs*n_per_rec
      if(.not.grid_written) irqlen=max(irqlen,nctrs*nqmax)
      if (.not.MA_Push_get(mt_dbl,irqlen,'rq',lrq,irq))
     &   call errquit('grid_den: cannot allocate rq',0, MA_ERR)
      if (.not.MA_Push_Get(mt_dbl,irqlen,'bp',lbp,ibp))
     &     call errquit('grid_quadv0: cannot allocate bp',4, MA_ERR)
c
c     compute variables for scratch space in function evaluation
c     assumes cartesian; spherical set would only require 2*MXang + 1
c     
      max_component = (nbf_ao_mxang+1)*(nbf_ao_mxang+2)/2 
c     
c     Allocate scratch space.
c     
      if (nbf_ao_mxang .eq. 0)then
         max_elem_ang_scr = max_component * 3
      elseif (nbf_ao_mxang .le. 3)then
         max_elem_ang_scr = max_component * 9
      else                      ! general case
         max_elem_ang_scr = max_component * 28
      endif 
c
c     Allocate and read in DIM coordinates
      i3 = 3*DIM_nAtoms
      if(.not.ma_push_get(mt_dbl,i3,'dimqm xyz',l_dimxyz,k_dimxyz))
     $   call errquit('dimqm_main malloc k_xyz failed',1,MA_ERR)
c   Read in coords from RTDB
      if(.not. rtdb_get(rtdb,'dimqm:coords',mt_dbl,i3,
     $                  dbl_mb(k_dimxyz)))
     $   call errquit('dimqm:xyz rtdb get failed',1,RTDB_ERR)
c
c     Allocate and read in induced dipoles      
      if(.not.ma_push_get(mt_dbl, i3*nmat, 'dimqm:muind',
     $                     l_muind, k_muind))
     $  call errquit('quadv0 malloc muind failed',2112, MA_ERR)

      call ga_init_fence()
      call ga_get(g_DIMmuind, 1, DIM_nOrder, 1, 1, dbl_mb(k_muind), 1)
      call ga_fence()

      call ga_sync()
c      
c     Grid quadrature
c
         if (iAOacc .le. 0) then
          if (.not. rtdb_get(rtdb, 'dft:iAOacc', mt_int, 1, iAOacc))then
             iAOacc=-nint(log(e_conv))
          else
             iAOacc=max(iAOacc,-nint(log(e_conv)))
          endif
         end if
c
         call dim_grid_gs_a(ncoef_max,nctrs, 
     &        int_mb(iniz), int_mb(ifin), 
     &        dbl_mb(irho),  
     &        dbl_mb(iqwght), dbl_mb(iqxyz), dbl_mb(ixyz), 
     &        dbl_mb(iexpo), dbl_mb(icoef), 
     &        dbl_mb(iBmat), dbl_mb(iFmat), 
     &        dbl_mb(iPmat), dbl_mb(irq), 
     &        int_mb(icetobfr),
     ,        int_mb(iibf_ao), dbl_mb(irqbuf),
     ,        dbl_mb(irchi_atom), g_vdim,
     ,        nmat, do_scf, do_cpks_l,
     ,        l3d, ipm, imag, dbl_mb(k_dimxyz),
     $        dbl_mb(k_muind), rtdb)
c     
      
      if (.not.ma_chop_stack(lqxyz))
     &   call errquit('grid_den: cannot pop stack',0, MA_ERR)
      call ga_sync()
      if(do_scf.and.(.not.l3d)) then
          call ga_symmetrize(g_vdim)
      endif
c
      if (.not.MA_Free_Heap(lxyz))
     &   call errquit('grid_den: cannot free heap',0, MA_ERR)
c
 123  format(
     &   1x,'Grid integrated density: ',f20.12,/
     &   1x,'Requested integration accuracy: ',e10.2)
c
C      if(oprint) then
C        write(luout,*) "dim_grid_quadv0 end"
C      end if

      return
      end
      Subroutine dim_grid_gs_a(ncoef_max, natoms, iniz, ifin,
     $                            rho, qwght, qxyz, xyz, expo, ccoef,
     $                            Bmat, Fmat, Pmat, rq, cetobfr,
     $                            ibf_ao, rqbuf, rchi_atom, g_vdim,
     $                            nmat, do_scf, do_cpks_l, l3d,
     $                            ipm, imag, dimxyz, muind, rtdb)

      implicit none
#include "errquit.fh"
#include "cgridfile.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "util.fh"
#include "steric.fh"
#include "msgids.fh"
#include "dimqm.fh"
c
c   Input Variables
      double precision rqbuf(*)
      double precision rad
      integer g_vdim
      double precision qxyz(3,n_per_rec) 
      double precision qwght(n_per_rec)
      integer ncoef_max
      integer iniz(natoms)
      integer ifin(natoms)
      double precision rho(*)
      double precision xyz(3,natoms)
      double precision expo(nbf_ao_mxprim)
      double precision Bmat(nbf_ao_mxnbf_ce*n_per_rec)
      double precision Pmat(*)
      double precision Fmat(*)
      double precision rq(n_per_rec, natoms)
      integer cetobfr(2,natoms)
      integer ibf_ao(nbf_ao)
      double precision rchi_atom(natoms)
      integer nmat
      logical do_scf
      logical do_cpks_l
      double precision ccoef(ncoef_max)
      integer natoms
      
      integer iqsh
      integer nqpts, ictr_buf
      integer ncube,istep,ntot_cube,ncontrset
      integer lbas_cent_info, ibas_cent_info,
     &        lbas_cset_info, ibas_cset_info,
     &        ldocset, idocset,i_iscratch,l_iscratch
      integer ncontrsetx,lbas_cent_xinfo, ibas_cent_xinfo,
     &        lbas_cset_xinfo, ibas_cset_xinfo,
     .     ldocsetx, idocsetx

      logical grid_file_rewind
      external grid_file_rewind
      logical xc_chkgrad, xc_chktau, kske
      external xc_chkgrad, xc_chktau
      logical l3d
      integer ipm
      integer imag
      double precision dimxyz(3, DIM_nAtoms)
      double precision muind(3, DIM_nAtoms)
      integer rtdb, id
      logical ldebug
c
      ldebug = .true.
      id = ga_nodeid()
c      if(ldebug.and.id.eq.0) then
c        write(luout,*) "dim_grid_quadv0a start"
c      end if
      
      if(n_rec_in_file.eq.0) goto 100
c
c     rewind grid pts file
c
      if (.not. grid_file_rewind())
     $   call errquit('grid_quadv0a: rewinding gridpts?', 0,
     &       UNKNOWN_ERR)
      if (.not.bas_numcont(AO_bas_han, ncontrset))
     &     call errquit('Exiting in grid_quadv0a',0, BASIS_ERR)
c
c     Allocate and create info for new basis function evaluator
c
      if (.not.MA_Push_Get(mt_int, 3*natoms, 'bas_cent_info',
     &     lbas_cent_info, ibas_cent_info))
     &     call errquit('grid_quadv0a: cannot allocate bas_cent_info',0,
     &       MA_ERR)
      if (.not.MA_Push_Get(mt_int, 6*ncontrset, 'bas_cset_info',
     &     lbas_cset_info, ibas_cset_info))
     &     call errquit('grid_quadv0a: cannot allocate bas_cset_info',0,
     &       MA_ERR)

      call xc_make_basis_info(AO_bas_han, int_mb(ibas_cent_info),
     &     int_mb(ibas_cset_info), natoms)

      if (.not.MA_Push_Get(mt_log, ncontrset, 'docset',
     &     ldocset, idocset))
     &     call errquit('grid_quadv0a: cannot allocate docset',
     .     ncontrset, MA_ERR)
c
      ntot_cube=0
c      write(luout,*) "n_recs_in_file:", n_rec_in_file
      do 200 iqsh = 1, n_rec_in_file
c      write(luout,*) "Top of loop", iqsh
c
c     Define the current range of radial shells and integration center.
c     
         call grid_file_read(n_per_rec, nqpts, ictr_buf, 
     &        rad,rqbuf,nsubb)
         if(nqpts.gt.buffer_size) call
     '        errquit(' buffersize exceed by qpts ',nqpts, DISK_ERR)
        if(nqpts.eq.0) goto 200
        istep=0
c        write(luout,*) iqsh, "nsubb:", nsubb
        do  ncube=1,nsubb
c     
c        temp put buf into currently used arrays qxyz and qwght
c     
         call grid_repack(rqbuf, qxyz, qwght, nqpts, rad,istep)
         if(nqpts.ne.0) then
            call dim_grid_gs_b(
     $              ictr_buf, qxyz, qwght, nqpts, rad,
     $              ncoef_max, natoms, iniz, ifin, rho,
     $              xyz, expo, Bmat, Fmat, Pmat, rq,
     $              cetobfr, ibf_ao, int_mb(ibas_cent_info),
     $              int_mb(ibas_cset_info), log_mb(idocset), rchi_atom,
     $              g_vdim, nmat, do_scf, do_cpks_l, l3d, ipm, imag,
     $              dimxyz, muind, rtdb)
c
             ntot_cube=ntot_cube+1
c
c
         endif
      enddo
c      write(luout,*) "Bottom of loop", iqsh
 200  continue
      call ga_igop(Msg_Excrho, ntot_cube , 1, '+')
c
      if (.not.ma_chop_stack(lbas_cent_info))
     &     call errquit('grid_quadv0a: pop stack failed.',1, MA_ERR)
 100  continue
c
c      if(ldebug) then
c        write(luout,*) "dim_grid_quadv0a end"
c      end if
      call ga_sync
      return
      end
      subroutine dim_grid_gs_b(ictr_buf, qxyz, qwght, nqpts, rad,
     $                            ncoef_max, natoms, iniz, ifin, rho,
     $                            xyz, expo, Bmat, Fmat, Pmat, rq,
     $                            cetobfr, ibf_ao, bas_cent_info,
     $                            bas_cset_info, docset, rchi_atom,
     $                            g_vdim, nmat, do_scf, do_cpks_l, l3d,
     $                            ipm, imag, dimxyz, muind, rtdb)
c
c$Id: dim_grid_quadv0b.F 19821 2010-12-14 07:46:49Z d3y133 $
c
      implicit none
#include "errquit.fh"
#include "cgridfile.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "grid_sic.fh"
#include "dftps.fh"
#include "geom.fh"
#include "dimqm.fh"
c
c   Input Variables
      integer ictr_buf ! ctr of grid
      double precision qxyz(3, n_per_rec) ! Quadrature point coordinates
      double precision qwght(n_per_rec)   ! Quadrature point weights
      integer nqpts                       ! Number of quadrature points
      double precision rad
      integer ncoef_max
      integer natoms    ! Number of atoms
      integer iniz(natoms)
      integer ifin(natoms)
      double precision rho(*)
      double precision xyz(3, natoms)   ! Atom coordinates
      double precision expo(nbf_ao_mxprim)
      double precision Bmat(nbf_ao_mxnbf_ce*n_per_rec)
      double precision Fmat(*)
      double precision Pmat(*)
      double precision rq(n_per_rec, natoms)
      integer cetobfr(2, natoms)
      integer ibf_ao(nbf_ao)
      integer bas_cent_info(3,natoms)
      integer bas_cset_info(6,*)
      logical docset(*)
      double precision rchi_atom(natoms)
      integer g_vdim ! DIM potential for this block of the quadrature
      integer nmat ! Number of density matrices (number of perturbing directions)
      logical do_scf
      logical do_cpks_l
      logical l3d
      integer ipm
      integer imag
      double precision dimxyz(3, DIM_nAtoms) ! DIM Coordinates
      double precision muind(3, DIM_nAtoms) ! DIM induced dipoles
      integer rtdb
c      
c   Local Variables 
      integer nbf_ao_mxnbf_ce2, mbf_ao, npol, ncontrset, maxdbas
      integer lchi_ao, ichi_ao, ldelchi_ao, idelchi_ao
      integer ldmat,idmat,i, dim_grid_nbfm
      double precision dabsmax
      double precision acc_ao_gauss
      external dabsmax
      external dim_grid_nbfm
      integer iscf_rho,iscf_delrho
      integer iscf_tau,iscf_ttau
      logical grid_sicinit,grid_sicend
      external grid_sicinit,grid_sicend
      double precision acc_xc_gauss
      integer lemat,iemat,lfmat,ifmat,k_scr,l_scr, id
      logical do_2nd
      logical stat
      logical dprint
      logical ldebug
      
      ldebug = .true.
      id = ga_nodeid()
c
      npol = (ipol*(ipol+1))/2
      acc_ao_gauss = dble(iaoacc)
      nbf_ao_mxnbf_ce2 = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
      if(nqpts.eq.0) return
      dprint = ldebug
c      if(dprint) then
c        write(luout,*) "dim_grid_quadv0b start"
c      end if
       maxdbas = 0
       idelchi_ao = 0
c
c     Evaluate the AO basis set at each of the quad. points.
c     allocate arrays for exponents and contraction coefficients
c     within int_eval_gbsets
c     Evaluate the distances (squared) between each center and the points
c
      mbf_ao=nbf_ao
      if(natoms.gt.1) then
         call icopy(mbf_ao, 0,0, ibf_ao, 1)  

         mbf_ao=dim_grid_nbfm(AO_bas_han,  natoms, 
     &        ictr_buf,rad,xyz,qxyz,nqpts,
     &        ibf_ao, docset, iniz, ifin, expo, 
     &        minexp,ldiff,acc_ao_gauss,iatype_pt_chg)
         if (mbf_ao.eq.0) return
      else
        if (.not.bas_numcont(AO_bas_han, ncontrset))
     &     call errquit('Exiting in grid_quadv0b',0, BASIS_ERR)
        iniz(1)=1
        ifin(1)=nbf_ao
        do i=1,ncontrset
          docset(i)=.true.
        enddo
        do i=1,nbf_ao
          ibf_ao(i)=i
        enddo
      endif

      if (.not.MA_Push_Get(mt_dbl, nqpts*mbf_ao, 'chi_ao', 
     &     lchi_ao, ichi_ao))
     &     call errquit('grid_quadv0b: cannot allocate chi_ao',0,
     &       MA_ERR)

      call qdist(rchi_atom, rq, qxyz, xyz, nqpts, natoms)

      call xc_eval_basis(ao_bas_han, maxdbas, dbl_mb(ichi_ao),
     &     0d0, 0d0, 0d0, rq, qxyz, xyz, nqpts, natoms,
     &     iniz, docset, bas_cent_info, bas_cset_info)
c
c     get reduced Xi(r) and dXi(r) over atoms
c
      call util_rmsatbf(nqpts, natoms,iniz,ifin,
     ,     dbl_mb(ichi_ao),rchi_atom)

c ---------------------
      if(nqpts.eq.0) goto 2010

c
c     Evaluate the DIM potential
c     Set up pointer to the SCF density for the CPKS LHS case.
c
c      call dim_eval_fnl(rho, nqpts, qwght,
c     &     dimxyz, qxyz, 
c     $     nmat, ipm, imag, muind,nDIM, rtdb)
      if (do_cpks_l) then
            call dim_eval_fnl(rho, nqpts, qwght, dimxyz,
     &            qxyz, nmat, ipm, imag, muind, DIM_nAtoms, rtdb)
      else
            call dim_eval_fnl_gs(rho, nqpts, qwght,
     &            dimxyz, qxyz, 
     $            muind,DIM_nAtoms, rtdb)
      endif
c       rho(1:nqpts*nmat) = 0.0d0
c
c   Transform the potential from over the quadrature points to the shape of the Fock Matrix
      if(do_scf) then
c
        call dim_tabcd(0, l3d, Fmat, Pmat, rho, Bmat,
     $                   dbl_mb(ichi_ao), 1, nqpts, mbf_ao,
     $                   nbf_ao_mxnbf_ce, nbf_ao_mxnbf_ce2,
     $                   AO_bas_han, natoms, iniz, ifin, g_vdim,
     $                   ibf_ao, cetobfr)
c
      else if (do_cpks_l) then
c
        call dim_tabcd(1, l3d, Fmat, Pmat, rho, Bmat,
     $                   dbl_mb(ichi_ao), nmat, nqpts, mbf_ao,
     $                   nbf_ao_mxnbf_ce, nbf_ao_mxnbf_ce2,
     $                   AO_bas_han, natoms, iniz, ifin, g_vdim,
     $                   ibf_ao, cetobfr)
c
      endif
c         write(luout,*) "g_dim AFTER dim_tabcd"
c         call ga_print(g_vdim)
c
c      endif

 2010 continue

      if (sic_orb_index.eq.1) then
         if(.not.grid_sicend(l_vect1,ldelrho_sig))
     ,        call errquit(' grid_quadv0b: grid_sicend failed',0,
     &       CALC_ERR)
      endif

      if (.not.ma_pop_stack(lchi_ao))
     &     call errquit('grid_quadv0b: cannot pop stack', 3, MA_ERR)

c      if(dprint) then
c        write(luout,*) "end dim_grid_quadv0b"
c      end if
      return
      end subroutine dim_grid_gs_b