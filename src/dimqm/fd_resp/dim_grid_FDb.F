      subroutine dim_grid_FDb(ictr_buf, qxyz, qwght, nqpts, rad,
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
      integer ipm, ii, jj
      integer imag
      double precision dimxyz(3, DIM_nAtoms) ! DIM Coordinates
      double precision muind(3, DIM_nAtoms, 3) ! DIM Coordinates
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

      if (DIM_ldebug .and. (id.eq.0)) 
     $ write(luout, *)'starting dim_grid_quadv0b'
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
      call dim_eval_fnl_FD(rho, nqpts, qwght, dimxyz,
     &            qxyz, nmat, muind, rtdb)
c       rho(1:nqpts*nmat) = 0.0d0
c
c      if (id.eq.0) then
c            do i = 1, nqpts
c                  write(luout,*)'ptnl: ', rho(i)
c            enddo
c      endif

c   Transform the potential from over the quadrature points to the shape of the Fock Matrix
c
      call dim_tabcd_FD(1, l3d, Fmat, Pmat, rho, Bmat,
     $                   dbl_mb(ichi_ao), nmat, nqpts, mbf_ao,
     $                   nbf_ao_mxnbf_ce, nbf_ao_mxnbf_ce2,
     $                   AO_bas_han, natoms, iniz, ifin, g_vdim,
     $                   ibf_ao, cetobfr)
c
c         write(luout,*) "g_dim AFTER dim_tabcd"
c         call ga_print(g_vxc)
c
c      endif

      !call ga_print(g_vdim)
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
      end subroutine dim_grid_FDb