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
c  NUFFT 1.3 release notes:
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
        subroutine nufft2d_kb_param(ms,eps,nf1,beta,nspread)
        implicit real *8 (a-h,o-z)
c        
        done=1
        pi=4*atan(done)
c
        rat = 2.5d0
        if( eps .le. 1d-4 ) nspread = 3
        if( eps .le. 1d-8 ) nspread = 5
        if( eps .le. 1d-12) nspread = 7
        if( eps .le. 1d-16) nspread = 9
c
        nf1 = rat*ms
        if (2*nspread.gt.nf1) then
        nf1 = next235(2d0*nspread) 
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
        subroutine nufft2d_kb_fm_init(nf1,beta,nspread,ms,fm)
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
        subroutine nufft2d_kb_xc_init(nf1,beta,nspread,nj,xj,xc)
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
      subroutine nufft2dkb0(nj,eps,ms,mt,lused)
      implicit real *8 (a-h,o-z)

      call nufft2d_kb_param(ms,eps,nf1,beta,nspread)
      call nufft2d_kb_param(mt,eps,nf2,beta,nspread)
        
      lused = 4 + (ms/2+1) + (mt/2+1) + 2*(2*nspread+1)*nj

      return
      end
c
c
c
c
      subroutine nufft2dkbi(nj,xj,yj,eps,ms,mt,ier,w,lw,lused)
      implicit real *8 (a-h,o-z)
      real*8 w(lw),xj(nj),yj(nj)

      ier=0

      if( lw .lt. 4) then
      ier=1
      return
      endif

      call nufft2d_kb_param(ms,eps,nf1,beta,nspread)
      call nufft2d_kb_param(mt,eps,nf2,beta,nspread)
        
      w(1)=nf1
      w(2)=nf2
      w(3)=beta
      w(4)=nspread
      iw1=1+4
      iw2=iw1+(ms/2+1)
      iw3=iw2+(mt/2+1)
      iw4=iw3+(2*nspread+1)*nj
      iw5=iw4+(2*nspread+1)*nj
      lused = 4 + (ms/2+1) + (mt/2+1) + 2*(2*nspread+1)*nj

      if( lused .gt. lw ) then
      ier=1
      return
      endif

      call nufft2d_kb_fm_init(nf1,beta,nspread,ms,w(iw1))
      call nufft2d_kb_fm_init(nf2,beta,nspread,mt,w(iw2))
      call nufft2d_kb_xc_init(nf1,beta,nspread,nj,xj,w(iw3))
      call nufft2d_kb_xc_init(nf2,beta,nspread,nj,yj,w(iw4))

      return
      end
c
c
c
c
      subroutine nufft2d1f90kbw(nj,xj,yj,cj,iflag,eps,ms,mt,fk,ier,w)
      implicit real *8 (a-h,o-z)
      real*8 w(*)

      nf1=w(1)
      nf2=w(2)
      beta=w(3)
      nspread=w(4)
      iw1=5
      iw2=iw1+(ms/2+1)
      iw3=iw2+(mt/2+1)
      iw4=iw3+(2*nspread+1)*nj
      iw5=iw4+(2*nspread+1)*nj

      call nufft2d1f90_kb(nj,xj,yj,cj,iflag,eps,ms,mt,fk,ier,
     $     nf1,nf2,nspread,w(iw3),w(iw4),w(iw1),w(iw2))
      
      return
      end
c
c
c
c
      subroutine nufft2d2f90kbw(nj,xj,yj,cj,iflag,eps,ms,mt,fk,ier,w)
      implicit real *8 (a-h,o-z)
      real*8 w(*)

      nf1=w(1)
      nf2=w(2)
      beta=w(3)
      nspread=w(4)
      iw1=5
      iw2=iw1+(ms/2+1)
      iw3=iw2+(mt/2+1)
      iw4=iw3+(2*nspread+1)*nj
      iw5=iw4+(2*nspread+1)*nj

cc        t1=second()
      call nufft2d2f90_kb(nj,xj,yj,cj,iflag,eps,ms,mt,fk,ier,
     $     nf1,nf2,nspread,w(iw3),w(iw4),w(iw1),w(iw2))
cc        t2=second()
cc        write(*,*) 'time',t2-t1
      return
      end
c
c
c
c
**********************************************************************
      subroutine nufft2d1f90_kb(nj,xj,yj,cj, iflag,eps, ms,mt,fk,ier,
     $      nf1,nf2,nspread,xc,yc,fm1,fm2)
      implicit none
      integer ier,iflag,ii,istart,iw10,iw11,iw13,iw14,iw15,iw16,iwtot
      integer j,jb1,jb2,jb1u,jb1d,k1,k2
      integer ms,mt,next235,nf1,nf2,nj,nspread,nw1,nw2,nw3
      real*8 rat,t1,t2,diff1,diff2,r2lamb,hx,hy
      real*8 cross,cross1,pi,eps
      real*8 xj(nj),yj(nj)
      real *8 xc(-nspread:nspread,nj),yc(-nspread:nspread,nj)
      real *8 fm1(0:ms/2),fm2(0:mt/2)
      parameter (pi=3.141592653589793238462643383279502884197d0)
      complex*16 cz,ccj,zz
      complex*16 cj(nj),fk(-ms/2:(ms-1)/2,-mt/2:(mt-1)/2)
c ----------------------------------------------------------------------
      real*8, allocatable :: fw(:)
      complex *16, allocatable :: cw(:)
c ----------------------------------------------------------------------
c
c     if (iflag .ge. 0) then
c
c               1  nj
c     fk(k1,k2) = -- SUM cj(j) exp(+i (k1,k2) * (xj(j, yj(j)) )   
c               nj j=1                     
c                                   for -ms/2 <= k1 <= (ms-1)/2
c                                   for -mt/2 <= k1 <= (mt-1)/2
c
c    else 
c
c               1  nj
c     fk(k1,k2) = -- SUM cj(j) exp(-i (k1,k2) * (xj(j, yj(j)) )   
c               nj j=1                     
c                                   for -ms/2 <= k1 <= (ms-1)/2
c                                   for -mt/2 <= k1 <= (mt-1)/2
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
c     xj,yj  location of sources (real *8)
c     cj     strengths of sources (complex *16)
c     iflag  determines sign of FFT (see above)
c     eps    precision request  (between 1.0d-33 and 1.0d-1)
c               recomended value is 1d-15 for double precision calculations
c     ms     number of Fourier modes computed (-ms/2 to (ms-1)/2 )
c     mt     number of Fourier modes computed (-mt/2 to (mt-1)/2 )
c
c     OUTPUT:
c
c     fk     Fourier transform values (2D complex *16 array)
c     ier    error return code
c   
c            ier = 0  => normal execution.
c            ier = 1  => precision eps requested is out of range.
c
c     The type 1 NUFFT proceeds in three steps (see [GL]).
c
c     1) spread data to oversampled regular mesh using convolution with
c        a Gaussian 
c     2) compute FFT on uniform mesh
c     3) deconvolve each Fourier mode independently
c          (mutiplying by Fourier transform of Gaussian).
c
c ----------------------------------------------------------------------
c
c     The oversampled regular mesh is defined by 
c
c     nf1 = rat*ms  points, where rat is the oversampling ratio.
c     nf2 = rat*mt  points, where rat is the oversampling ratio.
c       
c     For simplicity, we set  
c
c         rat = 2 for eps > 1.0d-11
c         rat = 3 for eps <= 1.0d-11.
c
c     The Gaussian used for convolution is:
c
c        g(x) = exp(-x^2 / 4tau) 
c
c     It can be shown [DR] that the precision eps is achieved when
c
c     nspread = int(-log(eps)/(pi*(rat-1d0)/(rat-.5d0)) + .5d0)
c     and tau is chosen as
c
c     tau = pi*lambda/(ms**2)
c     lambda = nspread/(rat(rat-0.5)).
c
c     Note that the Fourier transform of g(x) is
c
c     G(s) = exp(-s^2 tau) = exp(-pi*lambda s^2/ms^2)
c
c
c ----------------------------------------------------------------------
c     Fast Gaussian gridding is based on the following observation.
c
c     Let hx = 2*pi/nf1. In gridding data onto a regular mesh with
c     spacing nf1, we shift the source point xj by pi so 
c     that it lies in [0,2*pi] to simplify the calculations.
c     Since we are viewing the function
c     as periodic, this has no effect on the result.
c    
c     For source (xj+pi), let kb*hx denote the closest grid point and
c     let  kx*hx be a regular grid point within the spreading
c     distance. We can write
c
c     (xj+pi) - kx*hx = kb*hx + diff*hx - kx*hx = diff*hx - (kx-kb)*hx
c
c     where diff = (xj+pi)/hx - kb.
c
c     Let t1 = hx*hx/(4 tau) = pi/(nf1*nf1)/lambda*ms*ms
c                            = pi/lambda/(rat*rat)
c
c     exp(-( (xj+pi) -kx*hx)**2 / 4 tau)
c         = exp(-pi/lamb/rat^2 *(diff - (kx-kb))**2)
c         = exp(-t1 *(diff - (kx-kb))**2)
c         = exp(-t1*diff**2) * exp(2*t1*diff)**k * exp(-t1*k**2)
c           where k = kx-kb.
c 
c     In 2D/3D, this decomposition is carried out in each dimension.
c************************************************************************
c     -------------------------------
c     Precision dependent parameters
c
c     rat is oversampling parameter
c     nspread is number of neighbors to which Gaussian gridding is
c     carried out.
c     -------------------------------
      ier = 0
ccc        write(*,*) nspread
c
c     lambda (described above) = nspread/(rat*(rat-0.5d0)) 
c     It is more convenient to define r2lamb = rat*rat*lambda
c
      hx = 2*pi/nf1
      hy = 2*pi/nf2

c     -------------------------------
c     Compute workspace size and allocate
c
c     nw1 for fine grid, nw2 for second index, 
c     nw3 for max length to hold 1D FFT data.
c     -------------------------------
      nw1 = nf1 * nf2             
      nw2 = nw1 + nf2            
      nw3 = nw2 + max(nf1,nf2) 
c
      iw10 = 2*nw3
      iw11 = iw10 + ms
      iw13 = iw11 + mt/2 + 1
      iw14 = iw13 + nspread+1
      iw15 = iw14 + nspread+1
      iw16 = iw15 + 4*nf1+15
      iwtot = iw16 + 4*nf2+15
      allocate( fw(0:iwtot-1) )

      allocate( cw(0:nf1*nf2-1 + nf2) )
c
c     ---------------------------------------------------------------
c     Precompute spreading constants and initialize fw
c     to hold one of the terms in fast Gaussian gridding.
c     ---------------------------------------------------------------
      call zffti(nf1,fw(iw15))
      call zffti(nf2,fw(iw16))

c     ---------------------------------------------------------------
c     Initialize fine grid data to zero.
c     ---------------------------------------------------------------
      do k1 = 0, 2*nf1*nf2-1
         fw(k1) = 0.0d0
      enddo
      do k1 = 0, nf1*nf2-1
         cw(k1) = dcmplx(0.0d0,0.0d0)
      enddo
c
c     ---------------------------------------------------------------
c     Loop over sources (1,...,nj)
c
c     1. find closest mesh point (with periodic wrapping if necessary)
c     2. spread source data onto nearest nspread**2 grid points
c        using fast Gaussian gridding.
c
c     The following is a little hard to read because it takes
c     advantage of fast gridding and optimized to minimize the 
c     the number of multiplies in the inner loops.
c
c ---------------------------------------------------------------
      do j = 1, nj
         ccj = cj(j)/dble(nj)
         jb1 = int((xj(j)+pi)/hx)
         diff1 = (xj(j)+pi)/hx - jb1
         jb1 = mod(jb1, nf1)
         if (jb1.lt.0) jb1=jb1+nf1
         jb2 = int((yj(j)+pi)/hy)
         diff2 = (yj(j)+pi)/hy - jb2
         jb2 = mod(jb2, nf2)
         if (jb2.lt.0) jb2=jb2+nf2

         jb1d = min(nspread-1, jb1)
         jb1u = min(nspread, nf1-jb1-1)
         do k2 = -nspread+1, nspread
            ii = jb2+k2
            if (ii.lt.0) then
               ii = ii + nf2
            elseif (ii.ge.nf2) then
               ii = ii - nf2
            endif
            ii = jb1 + ii*nf1 
            cz = yc(k2,j)*ccj
            do k1 = -nspread+1, -jb1d-1
               istart = (k1+nf1+ii)
               cw(istart) = cw(istart) + xc(k1,j)*cz
            enddo
            do k1 = -jb1d, jb1u
               istart = (k1+ii)
               cw(istart) = cw(istart) + xc(k1,j)*cz
            enddo
            do k1 = jb1u+1, nspread
               istart = (k1-nf1+ii)
               cw(istart) = cw(istart) + xc(k1,j)*cz
            enddo
         enddo
      enddo
c
c     ---------------------------------------------------------------
c     Compute 2D FFT and carry out deconvolution.
c 
c     There is a factor of (-1)**k1 needed to account for the 
c     FFT phase shift.
c     ---------------------------------------------------------------
c
      do k2 = 0, nf1*nf2-1, nf1 
         if (iflag .ge. 0) then
            call zfftb(nf1,cw(k2),fw(iw15))
         else
            call zfftf(nf1,cw(k2),fw(iw15))
         endif
      enddo
c
      do k1 = -ms/2, (ms-1)/2
         ii = k1
         if (k1.lt.0) ii = nf1+k1
         do k2 = 0, nf2-1
            cw((nw1+k2)) = cw((ii + k2*nf1))
         enddo
         if (iflag .ge. 0) then
           call zfftb(nf2,cw(nw1),fw(iw16))
         else
           call zfftf(nf2,cw(nw1),fw(iw16))
         endif
c
         cross = fm1(abs(k1))
	 zz = cw(nw1)
         fk(k1, 0) = (cross*fm2(0))*zz
         do k2 = 1, (mt-1)/2
	    zz = cw(nw1+k2)
            fk(k1,k2) = (cross*fm2(k2))*zz
	    zz = cw(nw1+nf2-k2)
            fk(k1,-k2) = (cross*fm2(k2))*zz
         enddo
         if (mt/2*2.eq.mt) then 
            zz = cw(nw1+nf2-mt/2)
            fk(k1,-mt/2) = (cross*fm2(mt/2))*zz
         endif
      enddo
      deallocate(fw)
      deallocate(cw)
      return
      end
c
c
c
c
c
************************************************************************
      subroutine nufft2d2f90_kb(nj,xj,yj,cj, iflag,eps, ms,mt,fk,ier,
     $      nf1,nf2,nspread,xc,yc,fm1,fm2)
      implicit none
      integer iflag,ii,istart,ier,iw10,iw11,iw13,iw14,iw15,iw16,iwtot
      integer j,jb1,jb2,jb1u,jb1d,k1,k2
      integer ms,mt,next235,nf1,nf2,nj,nspread,nw1,nw2,nw3
      real*8 cross,cross1,diff1,diff2,eps,hx,hy,pi,rat,r2lamb,t1,t2
      real*8 xj(nj),yj(nj)
      real *8 xc(-nspread:nspread,nj),yc(-nspread:nspread,nj)
      real *8 fm1(0:ms/2),fm2(0:mt/2)
      parameter (pi=3.141592653589793238462643383279502884197d0)
      complex*16 cj(nj),fk(-ms/2:(ms-1)/2,-mt/2:(mt-1)/2),zz,cz
c ----------------------------------------------------------------------
      real*8, allocatable :: fw(:)
      complex *16, allocatable :: cw(:)
c ----------------------------------------------------------------------
c     if (iflag .ge. 0) then
c
c         (ms-1)/2 (mt-1)/2
c    cj(j) = SUM    SUM    fk(k1,k2) exp(+i (k1,k2) * (xj(j),yj(j)) )   
c         k1=-ms/2 k2=-mt/2             
c                                           for j = 1,...,nj
c     else
c
c         (ms-1)/2 (mt-1)/2
c    cj(j) = SUM    SUM    fk(k1,k2) exp(-i (k1,k2) * (xj(j),yj(j)) )   
c         k1=-ms/2 k2=-mt/2             
c                                           for j = 1,...,nj
c
c ----------------------------------------------------------------------
c     INPUT:
c
c     nj     number of output values   (integer)
c     xj,yj  location of output values (real *8 arrays)
c     iflag  determines sign of FFT (see above)
c     eps    precision request  (between 1.0d-33 and 1.0d-1)
c               recomended value is 1d-15 for double precision calculations
c     ms     number of Fourier modes given (1st index)  [-ms/2:(ms-1)/2]
c     mt     number of Fourier modes given (2nd index)  [-mt/2:(mt-1)/2]
c     fk     Fourier coefficient values (complex *16 2D array)
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
c     3) spread data to regular mesh using Gaussian
c
c
c     See subroutine nufft2d1f90
c     for more comments on fast gridding and parameter selection.
c
************************************************************************
c     -------------------------------
c     Precision dependent parameters
c
c     rat is oversampling parameter
c     nspread is number of neighbors to which Gaussian gridding is
c     carried out.
c     -------------------------------
      ier = 0
ccc        write(*,*) nspread
c
c     lambda (described above) = nspread/(rat*(rat-0.5d0)) 
c     It is more convenient to define r2lamb = rat*rat*lambda
c
c     -------------------------------

      hx = 2*pi/nf1
      hy = 2*pi/nf2
c
c     -------------------------------
c     Compute workspace size and allocate
c
c     nw1 for fine grid, nw2 for second index, 
c     nw3 for max length to hold 1D FFT data.
c     -------------------------------
      nw1 = nf1 * nf2
      nw2 = nw1 + nf2
      nw3 = nw2 + max(nf1,nf2)
c
      iw10 = 2*nw3
      iw11 = iw10 + ms
      iw13 = iw11 + mt/2 + 1
      iw14 = iw13 + nspread+1
      iw15 = iw14 + nspread+1
      iw16 = iw15 + 4*nf1+15
      iwtot = iw16 + 4*nf2+15
      allocate( fw(0:iwtot-1) )
      allocate( cw(0:nf1*nf2-1+nf2) )
c
c     ---------------------------------------------------------------
c     Precompute spreading constants and initialize ffts
c     ---------------------------------------------------------------

      call zffti(nf1,fw(iw15))
      call zffti(nf2,fw(iw16))
c
c
c     ---------------------------------------------------------------
c     Deconvolve and compute inverse 2D FFT
c     (A factor of (-1)**k is needed to shift phase.
c     ---------------------------------------------------------------
      do k1 = -ms/2, (ms-1)/2
         cross = fm1(abs(k1))
         zz  = (cross*fm2(0))*fk(k1,0)
         cw(nw1) = (zz)
         do k2 = 1, (mt-1)/2
            zz = (cross*fm2(k2))*fk(k1,k2)
            cw(nw1+k2) = (zz)
            zz = (cross*fm2(k2))*fk(k1,-k2)
            cw((nw1+nf2-k2)) = (zz)
         enddo
         if (mt/2*2.eq.mt) then
            zz = (cross*fm2(mt/2))*fk(k1,-mt/2)
            cw((nw1+nf2-mt/2)) = (zz)
	 endif
         do k2 = (mt+1)/2, nf2-mt/2-1
            cw(nw1+k2) = 0d0
         enddo
         if (iflag .ge. 0) then
            call zfftb(nf2,cw(nw1),fw(iw16))
         else
            call zfftf(nf2,cw(nw1),fw(iw16))
         endif
c
         ii = k1
         if (ii.lt.0) ii = nf1+k1
         do k2 = 0, nf2-1
            cw(ii+k2*nf1) = cw(nw1+k2)
         enddo
      enddo
c
      do k2 = 0, nf1*nf2-1, nf1
         do k1 = (ms+1)/2, nf1-ms/2-1
            cw((k1+k2)) = 0.0d0
         enddo
         if (iflag .ge. 0) then
            call zfftb(nf1,cw(k2),fw(iw15))
         else
            call zfftf(nf1,cw(k2),fw(iw15))
         endif
      enddo
c
c     ---------------------------------------------------------------
c     Loop over target points (1,...,nj)
c
c     1. find closest mesh point (with periodic wrapping if needed)
c     2. spread source data onto nearest nspread grid points
c        using fast Gaussian gridding.
c
c     The following is a little hard to read because it takes
c     advantage of fast gridding and optimized to minimize the 
c     the number of multiplies in the inner loops.
c     ---------------------------------------------------------------

      do j = 1, nj
         cj(j) = dcmplx(0d0,0d0)
         jb1 = int((xj(j)+pi)/hx)
         diff1 = (xj(j)+pi)/hx - jb1
         jb1 = mod(jb1, nf1)
         if (jb1.lt.0) jb1=jb1+nf1
         jb2 = int((yj(j)+pi)/hy)
         diff2 = (yj(j)+pi)/hy - jb2
         jb2 = mod(jb2, nf2)
         if (jb2.lt.0) jb2=jb2+nf2

         jb1d = min(nspread-1, jb1)
         jb1u = min(nspread, nf1-jb1-1)
         do k2 = -nspread+1, nspread
            ii = jb2+k2
            if (ii.lt.0) then
               ii = ii + nf2
            elseif (ii.ge.nf2) then
               ii = ii - nf2
            endif
            ii = jb1 + ii*nf1
            cz = dcmplx(0d0, 0d0)
            do k1 = -nspread+1, -jb1d-1
               zz = cw(k1+nf1+ii)
               cz = cz + xc(k1,j)*zz
            enddo
            do k1 = -jb1d, jb1u
               zz = cw(k1+ii)
               cz = cz + xc(k1,j)*zz
            enddo
            do k1 = jb1u+1, nspread
               zz = cw(k1-nf1+ii)
               cz = cz + xc(k1,j)*zz
            enddo
            cj(j) = cj(j) + yc(k2,j)*cz
         enddo
      enddo
      deallocate(fw)
      deallocate(cw)
      return
      end
c
c
c
c
c
