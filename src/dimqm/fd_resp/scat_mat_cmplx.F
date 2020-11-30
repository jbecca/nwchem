      subroutine scat_mat_cmplx(TT,T,Ni,Nj,mi,mj,ifirst,jfirst,ibf,jbf)
      implicit none
#include "stdio.fh"
#include "dimqm_constants.fh"
      integer mi, mj, Ni, Nj
      integer ifirst, jfirst
      integer ibf(mi), jbf(mj)
      double complex T(mi,mj), TT(Ni,Nj)
      integer i, j, jj
c       scatter into matrix
      call zfill(Ni*Nj, ZERO_C, TT, 1)
      do j = 1,mj
        jj=jbf(j)-jfirst+1
        do i=1,mi
            TT((ibf(i)-ifirst+1),jj) = T(i,j)
        enddo
      enddo

      return
      end subroutine scat_mat_cmplx

      subroutine dft_3dacc_cmplx(g_fock,buf,imat,imat2,
     $               ifirst,ilast,jfirst,jlast,ldin)
      implicit none
#include "dimqm_constants.fh"
      integer g_fock
      integer imat,imat2
      integer ldin
      integer ifirst,ilast,jfirst,jlast
      double complex buf(*)
c
      integer lo(3) !array of starting indices for ga
      integer hi(3) !array of ending indices for ga
      integer ld(2) !array specifying leading dimensions for buffer
      ld(1) = imat2-imat+1
      ld(2) = ldin
      lo(1) = imat
      hi(1) = imat2
      lo(2) = ifirst
      hi(2) = ilast
      lo(3) = jfirst
      hi(3) = jlast
      call nga_acc(g_fock, lo, hi, buf, ld, ONE_C)
      end subroutine dft_3dacc_cmplx

      subroutine transp_mat_cmplx(TT,T,Ni,Nj)
      implicit none
      integer Ni, Nj
      double complex T(Nj,Ni), TT(Ni,Nj)
      integer jj,ii 
      do jj=1, Ni
        do ii=1, Nj
            T(ii,jj) = TT(jj,ii)
        enddo
      enddo

      return
      end subroutine transp_mat_cmplx

