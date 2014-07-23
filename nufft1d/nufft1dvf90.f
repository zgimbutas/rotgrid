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
c
c       Vectorized nuFFT1d routines
c
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
      subroutine nufft1d1vf90(nvec,nj,xj,cj,iflag,eps,ms,fk,ier)
      implicit none
      integer ier,iflag,istart,iw1,iwtot,iwsav, kk,nvec
      integer j,jb1,jb1u,jb1d,k1,ms,next235,nf1,nj,nspread
      real*8 cross,cross1,diff1,eps,hx,pi,rat,r2lamb,t1,tau
      real*8 xc(-147:147),xj(nj)
      parameter (pi=3.141592653589793238462643383279502884197d0)
      complex*16 cj(nj,nvec),fk(-ms/2:(ms-1)/2,nvec),zz,ccj
c ----------------------------------------------------------------------
      complex*16, allocatable :: fw(:,:)
      real*8, allocatable :: fc(:)
      complex*16, allocatable :: wsave(:)
c ----------------------------------------------------------------------
c     if (iflag .ge. 0) then
c
c                 1  nj
c     fk(k1,l) = -- SUM cj(j,l) exp(+i k1 xj(j))  for -ms/2 <= k1 <= (ms-1)/2 
c                nj j=1                            
c
c     else
c
c                 1  nj
c     fk(k1,l) = -- SUM cj(j,l) exp(-i k1 xj(j))  for -ms/2 <= k1 <= (ms-1)/2 
c                nj j=1                            
c
c     l = 1,...,nvec
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
c     nvec     number of vectors to transform   (integer)
c     nj     number of sources   (integer)
c     xj     location of sources (real *8)
c
c            on interval [-pi,pi].
c
c     cj     strengths of sources (complex *16)
c     iflag  determines sign of FFT (see above)
c     eps    precision request  (between 1.0d-33 and 1.0d-1)
c               recomended value is 1d-15 for double precision calculations
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
c************************************************************************
c
c     Precision dependent parameters
c
c     rat is oversampling parameter
c     nspread is number of neighbors to which Gaussian gridding is
c     carried out.
c -------------------------------
      ier = 0
      if ((eps.lt.1d-33).or.(eps.gt.1d-1)) then
         ier = 1
         return
      endif
      if (eps.le.1d-11) then
         rat = 3.0d0
      else 
         rat = 2.0d0
      endif
      nspread = int(-log(eps)/(pi*(rat-1d0)/(rat-.5d0)) + .5d0)
      nf1 = rat*ms
      if (2*nspread.gt.nf1) then
         nf1 = next235(2d0*nspread) 
      endif 
c
c     lambda (described above) = nspread/(rat*(rat-0.5d0)) 
c     It is more convenient to define r2lamb = rat*rat*lambda
c
      r2lamb = rat*rat * nspread / (rat*(rat-.5d0))
      hx = 2*pi/nf1
c
c     -----------------------------------
c     Compute workspace size and allocate
c     -----------------------------------
      iw1 = 2*nf1
      iwtot = iw1
      allocate ( fw(0:iwtot,nvec) )
      allocate ( fc(nspread+1) )
      allocate ( wsave(4*nf1 + 15) )
c
c     ---------------------------------------------------------------
c     Precompute spreading constants and initialize fw
c     to hold one term needed for fast Gaussian gridding 
c     ---------------------------------------------------------------
      t1 = pi/r2lamb
      do k1 = 1, nspread
         fc(k1) = exp(-t1*k1**2)
      enddo
      call zffti(nf1,wsave)
c
c     ---------------------------------------------------------------
c     Initialize fine grid data to zero.
c     ---------------------------------------------------------------
      do kk=1,nvec
      do k1 = 0, nf1
         fw(k1,kk) = dcmplx(0d0,0d0)
      enddo
      enddo
c
c     ---------------------------------------------------------------
c     Loop over sources (1,...,nj)
c
c     1. find closest mesh point (with periodic wrapping if necessary)
c     2. spread source data onto nearest nspread grid points
c        using fast Gaussian gridding.
c
c     The following is a little hard to read because it takes
c     advantage of fast gridding and optimized to minimize the 
c     the number of multiplies in the inner loops.
c
c    ---------------------------------------------------------------
c
      do j = 1, nj
         jb1 = int((xj(j)+pi)/hx)
         diff1 = (xj(j)+pi)/hx - jb1
         jb1 = mod(jb1, nf1)
         if (jb1.lt.0) jb1=jb1+nf1
         xc(0) = exp(-t1*diff1**2)
         cross = xc(0)
         cross1 = exp(2d0*t1 * diff1)
         do k1 = 1, nspread
            cross = cross * cross1
            xc(k1) = fc(k1)*cross
         enddo
         cross = xc(0)
         cross1 = 1d0/cross1
         do k1 = 1, nspread-1
            cross = cross * cross1
            xc(-k1) = fc(k1)*cross
         enddo
c
         jb1d = min(nspread-1, jb1)
         jb1u = min(nspread, nf1-jb1-1)
         do kk=1,nvec
         ccj = cj(j,kk)/dble(nj)
         do k1 = -nspread+1, -jb1d-1
	    istart = jb1+k1+nf1
            zz=xc(k1)*ccj
            fw(istart,kk)=fw(istart,kk)+zz
         enddo
         do k1 = -jb1d, jb1u
	    istart = jb1+k1
            zz=xc(k1)*ccj
            fw(istart,kk)=fw(istart,kk)+zz
         enddo
         do k1 = jb1u+1, nspread
	    istart = jb1+k1-nf1
            zz=xc(k1)*ccj
            fw(istart,kk)=fw(istart,kk)+zz
         enddo
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
      do kk=1,nvec
      if (iflag .ge. 0) then
         call zfftb(nf1,fw(0,kk),wsave)
      else
         call zfftf(nf1,fw(0,kk),wsave)
      endif
      enddo
c
      tau = pi * r2lamb / dble(nf1)**2
      cross1 = 1d0/sqrt(r2lamb)
      do kk=1,nvec
      zz = fw(0,kk)
      fk(0,kk) = cross1*zz
      enddo
      do k1 = 1, (ms-1)/2
         cross1 = -cross1
         cross = cross1*exp(tau*dble(k1)**2)
         do kk=1,nvec
	 zz = fw(k1,kk)
         fk(k1,kk) = cross*zz
	 zz = fw(nf1-k1,kk)
         fk(-k1,kk) = cross*zz
         enddo
      enddo
      if (ms/2*2.eq.ms) then
         cross = -cross1*exp(tau*dble(ms/2)**2)
         do kk=1,nvec
         zz = fw(nf1-ms/2,kk)
         fk(-ms/2,kk) = cross*zz
         enddo
      endif
      deallocate(fw)
      deallocate(fc)
      deallocate(wsave)
      return
      end
c
c
c
c
c
************************************************************************
      subroutine nufft1d2vf90(nvec,nj,xj,cj, iflag,eps, ms,fk,ier)
      implicit none
      integer ier,iflag,iw1,iwsav,iwtot,j,jb1,jb1u,jb1d,k1, kk,nvec
      integer ms,next235,nf1,nj,nspread,nw
      real*8 cross,cross1,diff1,eps,hx,pi,rat,r2lamb,t1
      real*8 xj(nj),xc(-147:147)
      parameter (pi=3.141592653589793238462643383279502884197d0)
      complex*16 cj(nj,nvec), fk(-ms/2:(ms-1)/2,nvec)
      complex*16 zz
c ----------------------------------------------------------------------
      complex*16, allocatable :: fw(:,:)
      real*8, allocatable :: fc(:)
      complex*16, allocatable :: wsave(:)
c ----------------------------------------------------------------------
c     if (iflag .ge. 0) then
c
c                (ms-1)/2
c     cj(j,l) =    SUM      fk(k1,l) exp(+i k1 xj(j))  for j = 1,...,nj
c                k1= -ms/2                            
c
c     else
c
c                (ms-1)/2
c     cj(j,l) =    SUM      fk(k1,l) exp(-i k1 xj(j))  for j = 1,...,nj
c                k1= -ms/2                            
c
c     l = 1,...,nvec
c
c ----------------------------------------------------------------------
c     INPUT:
c
c     nvec     number of vectors to transform   (integer)
c     nj     number of output values   (integer)
c     xj     location of output values (real *8 array)
c     iflag  determines sign of FFT (see above)
c     eps    precision request  (between 1.0d-33 and 1.0d-1)
c               recomended value is 1d-15 for double precision calculations
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
c     3) spread data to regular mesh using Gaussian
c
c
c     See subroutine nufft1d1f90(nj,xj,cj,iflag,eps,ms,fk,ier)
c     for more comments on fast gridding and parameter selection.
c
************************************************************************
c
c     Precision dependent parameters
c
c     rat is oversampling parameter
c     nspread is number of neighbors to which Gaussian gridding is
c     carried out.
c     -------------------------------
c
c     lambda (described above) = nspread/(rat*(rat-0.5d0)) 
c     It is more convenient to define r2lamb = rat*rat*lambda
c
c     -------------------------------
      ier = 0
      if ((eps.lt.1d-33).or.(eps.gt.1d-1)) then
         ier = 1
         return
      endif
      if (eps.le.1d-11) then
         rat = 3.0d0
      else 
         rat = 2.0d0 
      endif
c
      nspread = int(-log(eps)/(pi*(rat-1d0)/(rat-.5d0)) + .5d0)
      nf1 = rat*ms
      if (2*nspread.gt.nf1) then
         nf1 = next235(2d0*nspread) 
      endif 
c
      r2lamb = rat*rat * nspread / (rat*(rat-.5d0))
      hx = 2*pi/nf1
c
c     -----------------------------------
c     Compute workspace size and allocate
c     -----------------------------------
      iw1 = 2*nf1
      iwtot = iw1
      allocate ( fw(0:iwtot,nvec) )
      allocate ( fc(nspread+1) )
      allocate ( wsave(4*nf1 + 15) )
c
c     ---------------------------------------------------------------
c     Precompute spreading constants and initialize fw
c     to hold one term needed for fast Gaussian gridding 
c     ---------------------------------------------------------------
      t1 = pi/r2lamb
      do k1 = 1, nspread
         fc(k1) = exp(-t1*k1**2)
      enddo
      call zffti(nf1,wsave)
c
c     ---------------------------------------------------------------
c     Deconvolve and compute inverse 1D FFT
c     (A factor of (-1)**k is needed to shift phase.)
c     ---------------------------------------------------------------
c
      t1 = pi * r2lamb / dble(nf1)**2
      cross1 = 1d0/sqrt(r2lamb)
      do kk=1,nvec
      zz = cross1*fk(0,kk)
      fw(0,kk) = zz
      enddo
      do k1 = 1, (ms-1)/2
         cross1 = -cross1
         cross = cross1*exp(t1*dble(k1)**2)
         do kk=1,nvec
         zz = cross*fk(k1,kk)
         fw(k1,kk) = zz
         zz = cross*fk(-k1,kk)
         fw(nf1-k1,kk) = zz
         enddo
      enddo
      cross = -cross1*exp(t1*dble(ms/2)**2)
      if (ms/2*2.eq.ms) then
         do kk=1,nvec
	 zz = cross*fk(-ms/2,kk)
         fw(nf1-ms/2,kk) = zz
         enddo
      endif
      do k1 = (ms+1)/2, nf1-ms/2-1
         do kk=1,nvec
         fw(k1,kk) = dcmplx(0d0, 0d0)
         enddo
      enddo
c
      do kk=1,nvec
      if (iflag .ge. 0) then
         call zfftb(nf1,fw(0,kk),wsave)
      else
         call zfftf(nf1,fw(0,kk),wsave)
      endif
      enddo
c
c     ---------------------------------------------------------------
c     Loop over target points (1,...,nj)
c
c       1. find closest mesh point (with periodic wrapping if needed)
c       2. get contributions from regular fine grid to target
c          locations using Gaussian convolution.
c     ---------------------------------------------------------------
      t1 = pi/r2lamb
      do j = 1, nj
         jb1 = int((xj(j)+pi)/hx)
         diff1 = (xj(j)+pi)/hx - jb1
         jb1 = mod(jb1, nf1)
         if (jb1.lt.0) jb1=jb1+nf1
         xc(0) = exp(-t1*diff1**2)
         cross = xc(0)
         cross1 = exp(2d0*t1 * diff1)
         do k1 = 1, nspread
            cross = cross * cross1
            xc(k1) = fc(k1)*cross
         enddo
         cross = xc(0)
         cross1 = 1d0/cross1
         do k1 = 1, nspread-1
            cross = cross * cross1
            xc(-k1) = fc(k1)*cross
         enddo
c
         jb1d = min(nspread-1, jb1)
         jb1u = min(nspread, nf1-jb1-1)
         do kk=1,nvec
         cj(j,kk) = dcmplx(0d0,0d0)
         do k1 = -nspread+1, -jb1d-1
            zz = fw(jb1+k1+nf1,kk)
            cj(j,kk) = cj(j,kk) + xc(k1)*zz
         enddo
         do k1 = -jb1d, jb1u
            zz = fw(jb1+k1,kk)
            cj(j,kk) = cj(j,kk) + xc(k1)*zz
         enddo
         do k1 = jb1u+1, nspread
            zz = fw(jb1+k1-nf1,kk)
            cj(j,kk) = cj(j,kk) + xc(k1)*zz
         enddo
         enddo
      enddo
      deallocate(fw)
      deallocate(fc)
      deallocate(wsave)
      return
      end
c
c
c
c
c
c
        subroutine nufft1d3vf90(nvec,nj,xj,cj,iflag,eps,nk,sk,fk,ier)
        implicit real *8 (a-h,o-z)
        real*8 xj(nj), sk(nk)
        complex*16 cj(nj,nvec), fk(nk,nvec)

cccC$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k)
        do k=1,nvec
        call nufft1d3f90(nj,xj,cj(1,k),iflag,eps,nk,sk,fk(1,k),ier)
        enddo
cccC$OMP END PARALLEL DO

        return
        end
