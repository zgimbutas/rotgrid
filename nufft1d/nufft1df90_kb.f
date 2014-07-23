cc Copyright (C) 2004-2009: Leslie Greengard and June-Yub Lee 
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
cc
c
c  NUFFT 1.5 release notes:
c
c  These codes are asymptotically fast (O(N log N)), but not optimized.
c
c  1) We initialize the FFT on every call.
c
c  2) We do not precompute the exponentials involved in "fast Gaussian
c  gridding".
c
c  3) We do not block structure the code so that irregularly placed points
c  are interpolated (gridded) in a cache-aware fashion.
c
c  4) We use the Netlib FFT library (www.netlib.org) 
c     rather than the state of the art FFTW package (www.fftw.org).
c
c  Different applications have different needs, and we have chosen
c  to provide the simplest code as a reasonable efficient template.
c
c**********************************************************************
c
c
c
c
c
        subroutine nufft1d_kb_param(ms,eps,nf1,beta,nspread)
        implicit real *8 (a-h,o-z)
c        
        done=1
        pi=4*atan(done)
c
        rat = 2.5d0
        if( eps .le. 1d-4 ) nspread = 3+1
        if( eps .le. 1d-8 ) nspread = 5+1
        if( eps .le. 1d-12) nspread = 7+2
        if( eps .le. 1d-16) nspread = 9+2
c
        nf1 = rat*ms
        if (2*nspread.gt.nf1) then
        nf1 = next235(2d0*nspread) 
        else
        nf1 = next235(1d0*nf1) 
        endif 
c
        hx = 2*pi/nf1
c
        beta = (2-1.0d0/rat)*pi
c
        return
        end
c
c
c
c
c
        subroutine nufft1d_kb_fm_init(nf1,beta,nspread,ms,fm)
        implicit real *8 (a-h,o-z)
        real *8 fm(0:ms/2)

        done=1
        pi=4*atan(done)

        call i0eva(beta*nspread,fc)

        k1 = 0
        cross1 = 1
        tc = nspread * sqrt(beta**2-(2*pi*k1/dble(nf1))**2)
        call i0eva(tc,cross)
        cross = cross / fc * pi * cross1
        cross = 1/cross
        fm(0) = cross 
        
        do k1 = 1, (ms-1)/2
        cross1 = -cross1
        tc = nspread * sqrt(beta**2-(2*pi*k1/dble(nf1))**2)
        call i0eva(tc,cross)        
        cross = cross / fc * pi * cross1
        cross = 1/cross
        fm(k1) = cross
        enddo

        if (ms/2*2.eq.ms) then
        k1=ms/2
        cross1 = -cross1
        tc = nspread * sqrt(beta**2-(2*pi*k1/dble(nf1))**2)
        call i0eva(tc,cross)        
        cross = cross / fc * pi * cross1
        cross = 1/cross
        fm(ms/2) = cross
        endif
c
        return
        end
c
c
c
c
c
        subroutine nufft1d_kb_xc_init(nf1,beta,nspread,nj,xj,xc)
        implicit real *8 (a-h,o-z)
        real *8 xj(nj)
        real *8 xc(-nspread:nspread,nj)

        done=1
        pi=4*atan(done)

        hx=2*pi/nf1

        call i0eva(beta*nspread,fc)

        do j = 1,nj

        jb1 = int((xj(j)+pi)/hx)
        diff1 = (xj(j)+pi)/hx - jb1
        jb1 = mod(jb1, nf1)
        if (jb1.lt.0) jb1=jb1+nf1

        do k1 = -nspread+1,nspread
        tc = sqrt(nspread**2-(-k1+diff1)**2)
c        
c        if( tc .ne. 0 ) xc(k1,j) = sinh(beta*tc)/tc
c        if( tc .ne. 0 ) xc(k1,j) = (exp(beta*tc)-exp(-beta*tc))/2/tc
c
        if( tc .ne. 0 ) then
        ec = exp(beta*tc)
        xc(k1,j) = (ec-1/ec)/2/tc
        endif

        if( tc .eq. 0 ) xc(k1,j) = 1

        xc(k1,j)=xc(k1,j)/fc         
        enddo
        xc(-nspread,j)=0
         
        enddo
c
        return
        end
c
c
c
c
c**********************************************************************
c
c
c
c
      subroutine nufft1d1f90kb0(nj,eps,ms,lused)
      implicit real *8 (a-h,o-z)

      call nufft1d_kb_param(ms,eps,nf1,beta,nspread)
        
      iw1=1+3
      iw2=iw1+(ms/2+1)
      iw3=iw2+(2*nspread+1)*nj
      lused = 3 + (ms/2+1) + (2*nspread+1)*nj

      return
      end
c
c
c
c
      subroutine nufft1dkbi(nj,xj,eps,ms,ier,w,lw,lused)
      implicit real *8 (a-h,o-z)
      real*8 w(lw),xj(nj)

      ier=0

      if( lw .lt. 3) then
      ier=1
      return
      endif

      call nufft1d_kb_param(ms,eps,nf1,beta,nspread)
        
      w(1)=nf1
      w(2)=beta
      w(3)=nspread
      iw1=1+3
      iw2=iw1+(ms/2+1)
      iw3=iw2+(2*nspread+1)*nj
      lused = 3 + (ms/2+1) + (2*nspread+1)*nj

      if( lused .gt. lw ) then
      ier=1
      return
      endif

      call nufft1d_kb_fm_init(nf1,beta,nspread,ms,w(iw1))
      call nufft1d_kb_xc_init(nf1,beta,nspread,nj,xj,w(iw2))

      return
      end
c
c
c
c
      subroutine nufft1df90kbw(nj,xj,cj,iflag,eps,ms,fk,ier,w)
      implicit real *8 (a-h,o-z)
      real*8 w(*)

      nf1=w(1)
      beta=w(2)
      nspread=w(3)
      iw1=4
      iw2=iw1+(ms/2+1)
      iw3=iw2+(2*nspread+1)*nj

      call nufft1d1f90_kb(nj,xj,cj,iflag,eps,ms,fk,ier,
     $     nf1,nspread,w(iw1),w(iw2))
      
      return
      end
c
c
c
c
      subroutine nufft1d2f90kbw(nj,xj,cj,iflag,eps,ms,fk,ier,w)
      implicit real *8 (a-h,o-z)
      real*8 w(*)

      nf1=w(1)
      beta=w(2)
      nspread=w(3)
      iw1=4
      iw2=iw1+(ms/2+1)
      iw3=iw2+(2*nspread+1)*nj

      call nufft1d2f90_kb(nj,xj,cj,iflag,eps,ms,fk,ier,
     $     nf1,nspread,w(iw1),w(iw2))
      
      return
      end
c
c
c
c
c**********************************************************************
      subroutine nufft1d1f90_kb(nj,xj,cj,iflag,eps,ms,fk,ier,
     $     nf1,nspread,xc,fm)
      implicit real *8 (a-h,o-z)
      integer ier,iflag,istart,iw1,iwtot,iwsav
      integer j,jb1,jb1u,jb1d,k1,ms,next235,nf1,nj,nspread
      real*8 cross,cross1,diff1,eps,hx,pi,rat,r2lamb,t1,tau
      real*8 xc(-nspread:nspread,nj),fm(0:ms/2),xj(nj),beta
      parameter (pi=3.141592653589793238462643383279502884197d0)
      complex*16 cj(nj),fk(-ms/2:(ms-1)/2),zz,ccj
c ----------------------------------------------------------------------
      real*8, allocatable :: fw(:)
c ----------------------------------------------------------------------
c     if (iflag .ge. 0) then
c
c               1  nj
c     fk(k1) = -- SUM cj(j) exp(+i k1 xj(j))  for -ms/2 <= k1 <= (ms-1)/2 
c              nj j=1                            
c
c     else
c
c               1  nj
c     fk(k1) = -- SUM cj(j) exp(-i k1 xj(j))  for -ms/2 <= k1 <= (ms-1)/2 
c              nj j=1                            
c
c     References:
c
c     [DR] Fast Fourier transforms for nonequispaced data,
c          A. Dutt and V. Rokhlin, SIAM J. Sci. Comput. 14, 
c          1368-1383, 1993.
c
c     [GL] Accelerating the Nonuniform Fast Fourier Transform,
c          L. Greengard and J.-Y. Lee, SIAM Review 46, 443-454 (2004).
c
c ----------------------------------------------------------------------
c     INPUT:
c
c     nj     number of sources   (integer)
c     xj     location of sources (real *8)
c
c            on interval [-pi,pi].
c
c     cj     strengths of sources (complex *16)
c     iflag  determines sign of FFT (see above)
c     eps    precision request  (between 1.0d-33 and 1.0d-1)
c     ms     number of Fourier modes computed (-ms/2 to (ms-1)/2 )
c
c     OUTPUT:
c
c     fk     Fourier transform values (complex *16)
c     ier    error return code
c   
c            ier = 0  => normal execution.
c            ier = 1  => precision eps requested is out of range.
c
c     The type 1 NUFFT proceeds in three steps (see [GL]).
c
c     1) spread data to oversampled regular mesh using convolution
c     2) compute FFT on uniform mesh
c     3) deconvolve each Fourier mode independently
c
c ----------------------------------------------------------------------
c
      hx = 2*pi/nf1
c
c     -----------------------------------
c     Compute workspace size and allocate
c     -----------------------------------
      iw1 = 2*nf1
      iwsav = iw1+nspread+1
      iwtot = iwsav+4*nf1+15
      allocate ( fw(0:iwtot) )
c
      call zffti(nf1,fw(iwsav))
c
c     ---------------------------------------------------------------
c     Initialize fine grid data to zero.
c     ---------------------------------------------------------------
      do k1 = 0, 2*nf1-1
         fw(k1) = 0
      enddo
c
c     ---------------------------------------------------------------
c     Loop over sources (1,...,nj)
c
c     1. find closest mesh point (with periodic wrapping if necessary)
c     2. spread source data onto nearest nspread grid points
c    ---------------------------------------------------------------
c
      do j = 1, nj
         ccj = cj(j)/dble(nj)

         jb1 = int((xj(j)+pi)/hx)
         diff1 = (xj(j)+pi)/hx - jb1
         jb1 = mod(jb1, nf1)
         if (jb1.lt.0) jb1=jb1+nf1
c
         jb1d = min(nspread-1, jb1)
         jb1u = min(nspread, nf1-jb1-1)
         do k1 = -nspread+1, -jb1d-1
	    istart = 2*(jb1+k1+nf1)
            zz=xc(k1,j)*ccj
            fw(istart)=fw(istart)+dreal(zz)
            fw(istart+1)=fw(istart+1)+dimag(zz)
         enddo
         do k1 = -jb1d, jb1u
	    istart = 2*(jb1+k1)
            zz=xc(k1,j)*ccj
            fw(istart)=fw(istart)+dreal(zz)
            fw(istart+1)=fw(istart+1)+dimag(zz)
         enddo
         do k1 = jb1u+1, nspread
	    istart = 2*(jb1+k1-nf1)
            zz=xc(k1,j)*ccj
            fw(istart)=fw(istart)+dreal(zz)
            fw(istart+1)=fw(istart+1)+dimag(zz)
         enddo
      enddo
c
c     ---------------------------------------------------------------
c     Compute 1D FFT and carry out deconvolution.
c 
c     There is a factor of (-1)**k1 needed to account for the 
c     FFT phase shift.
c     ---------------------------------------------------------------
c
c
      if (iflag .ge. 0) then
         call zfftb(nf1,fw(0),fw(iwsav))
      else
         call zfftf(nf1,fw(0),fw(iwsav))
      endif
c
c
      zz = dcmplx(fw(0),fw(1))
      fk(0) = fm(0) * zz

      do k1 = 1, (ms-1)/2
        zz = dcmplx(fw(2*k1),fw(2*k1+1))
        fk(k1) = fm(k1)*zz
        zz = dcmplx(fw(2*(nf1-k1)),fw(2*(nf1-k1)+1))
        fk(-k1) = fm(k1)*zz
      enddo
      if (ms/2*2.eq.ms) then
         zz = dcmplx(fw(2*nf1-ms),fw(2*nf1-ms+1))
         fk(-ms/2) = fm(ms/2)*zz
      endif

      deallocate(fw)

      return
      end
c
c
c
c
c
************************************************************************
      subroutine nufft1d2f90_kb(nj,xj,cj, iflag,eps, ms,fk,ier,
     $     nf1,nspread,xc,fm)
      implicit real *8 (a-h,o-z)
      integer ier,iflag,iw1,iwsav,iwtot,j,jb1,jb1u,jb1d,k1
      integer ms,next235,nf1,nj,nspread,nw
      real*8 cross,cross1,diff1,eps,hx,pi,rat,r2lamb,t1,beta
      real*8 xj(nj),xc(-nspread:nspread,nj),fm(0:ms/2)
      parameter (pi=3.141592653589793238462643383279502884197d0)
      complex*16 cj(nj), fk(-ms/2:(ms-1)/2)
      complex*16 zz
c ----------------------------------------------------------------------
      real*8, allocatable :: fw(:)
c ----------------------------------------------------------------------
c     if (iflag .ge. 0) then
c
c              (ms-1)/2
c     cj(j) =    SUM      fk(k1) exp(+i k1 xj(j))  for j = 1,...,nj
c              k1= -ms/2                            
c
c     else
c
c              (ms-1)/2
c     cj(j) =    SUM      fk(k1) exp(-i k1 xj(j))  for j = 1,...,nj
c              k1= -ms/2                            
c
c ----------------------------------------------------------------------
c     INPUT:
c
c     nj     number of output values   (integer)
c     xj     location of output values (real *8 array)
c     iflag  determines sign of FFT (see above)
c     eps    precision request  (between 1.0d-12 and 1.0d-1)
c     ms     number of Fourier modes given  [ -ms/2: (ms-1)/2 ]
c     fk     Fourier coefficient values (complex *16 array)
c
c     OUTPUT:
c
c     cj     output values (complex *16 array)
c     ier    error return code
c   
c            ier = 0  => normal execution.
c            ier = 1  => precision eps requested is out of range.
c
c
c     The type 2 algorithm proceeds in three steps (see [GL]).
c
c     1) deconvolve (amplify) each Fourier mode first 
c     2) compute inverse FFT on uniform fine grid
c     3) spread data to regular mesh
c
************************************************************************
c
      hx = 2*pi/nf1
c
c     -----------------------------------
c     Compute workspace size and allocate
c     -----------------------------------
      iw1 = 2*nf1
      iwsav = iw1 + nspread+1
      iwtot = iwsav + 4*nf1 + 15
      allocate ( fw(0:iwtot))
c
      call zffti(nf1,fw(iwsav))
c
c     ---------------------------------------------------------------
c     Deconvolve and compute inverse 1D FFT
c     (A factor of (-1)**k is needed to shift phase.)
c     ---------------------------------------------------------------
c
      zz = fm(0)*fk(0)
      fw(0) = dreal(zz)
      fw(1) = dimag(zz)
      do k1 = 1, (ms-1)/2
         zz = fm(k1)*fk(k1)
         fw(2*k1) = dreal(zz)
         fw(2*k1+1) = dimag(zz)
         zz = fm(k1)*fk(-k1)
         fw(2*(nf1-k1)) = dreal(zz)
         fw(2*(nf1-k1)+1) = dimag(zz)
      enddo
      if (ms/2*2.eq.ms) then
	 zz = fm(ms/2)*fk(-ms/2)
         fw(2*nf1-ms) = dreal(zz)
         fw(2*nf1-ms+1) = dimag(zz)
      endif
      do k1 = (ms+1)/2, nf1-ms/2-1
         fw(2*k1) = dcmplx(0d0, 0d0)
         fw(2*k1+1) = dcmplx(0d0, 0d0)
      enddo
c
c
      if (iflag .ge. 0) then
         call zfftb(nf1,fw(0),fw(iwsav))
      else
         call zfftf(nf1,fw(0),fw(iwsav))
      endif
c
c     ---------------------------------------------------------------
c     Loop over target points (1,...,nj)
c
c       1. find closest mesh point (with periodic wrapping if needed)
c       2. get contributions from regular fine grid to target locations
c     ---------------------------------------------------------------
      do j = 1, nj
         cj(j) = dcmplx(0d0,0d0)
         jb1 = int((xj(j)+pi)/hx)
         diff1 = (xj(j)+pi)/hx - jb1
         jb1 = mod(jb1, nf1)
         if (jb1.lt.0) jb1=jb1+nf1

         jb1d = min(nspread-1, jb1)
         jb1u = min(nspread, nf1-jb1-1)
         do k1 = -nspread+1, -jb1d-1
	    zz = dcmplx(fw(2*(jb1+k1+nf1)),fw(2*(jb1+k1+nf1)+1))
            cj(j) = cj(j) + xc(k1,j)*zz
         enddo
         do k1 = -jb1d, jb1u
	    zz = dcmplx(fw(2*(jb1+k1)),fw(2*(jb1+k1)+1))
            cj(j) = cj(j) + xc(k1,j)*zz
         enddo
         do k1 = jb1u+1, nspread
	    zz = dcmplx(fw(2*(jb1+k1-nf1)),fw(2*(jb1+k1-nf1)+1))
            cj(j) = cj(j) + xc(k1,j)*zz
         enddo
      enddo
c
      deallocate(fw)
      return
      end
c
c
c
c
c
c
