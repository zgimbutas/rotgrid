cc Copyright (C) 2009-2012: Leslie Greengard and Zydrunas Gimbutas
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
c
c
c    $Date$
c    $Revision$
c
c
c     WIGNER ROTATION OPERATOR VIA PROJECTION  (FORTRAN 77 and 90 VERSIONS).
c
c     Requires FFT and Associated Legendre Function Libraries.
c
c     User-callable f77 routine is rotmat_proj. 
c     User-callable f90 routine is rotmat_projf90.
c     The other routines are used internally.
c
c***********************************************************************
      subroutine rotmat_proj(beta,alpha,nterms,m1,m2,mpole,lmp,
     $     marray2,lmpn,w,lw,lused)
c***********************************************************************
c       Purpose:
c
c	Fast and stable algorithm for constructing the rotation operator about
c	the z-axis determined by angle alpha plus the y-axis determined
c	by angle beta. After rotation, the expansion pole is moved to
c	location (beta, alpha) in spherical coordinates (theta, phi).
c
c       The method is based on computing the induced potential and
c       its theta-derivative on the rotated equator
c       for each order (first index). The coefficients of  the rotated
c       expansion can then be obtained by FFT and projection.
c
c       There is some loss in speed over using recurrence relations 
c       but it is stable to all orders whereas the recurrence schemes 
c       are not.
c
c       Our definition of complex spherical harmonics is
c
c       Ynm(theta,phi)= sqrt( 2n+1) sqrt((n-m)!/(n+m)!) 
c                       Pnm(cos theta) e^(im phi), 
c       Yn,-m(theta,phi) = sqrt( 2n+1) sqrt((n-m)!/(n+m)!) 
c                       Pnm(cos theta) e^(-im phi),   for m >= 0.
c       
c       Note that we do not include the Condon-Shortley phase (-1)^m
c
C---------------------------------------------------------------------
c       INPUT:
c
c       beta:  the rotation angle about the y-axis.
c       alpha:  the rotation angle about the z-axis.
c       nterms: order of multipole expansion
c
C       lmpn    leading dim for marray2 (must exceed nterms)
c       w     :  work array 
c       lw    :  length of work array 
c
C---------------------------------------------------------------------
c       OUTPUT:
c
c       marray2  coefficients of rotated expansion.
c       lused    amount of workspace used.
c       ier      error return flag
c                0 successful execution
c                1 insufficient memory
c
C---------------------------------------------------------------------
c
c
c
      implicit real *8 (a-h,o-z)
      integer nquad
      real *8 w(lw)
      complex *16 marray2(0:lmpn,-lmpn:lmpn,-lmpn:lmpn)
c
c      nquad = 2*nterms+2
c        write(*,*) nquad
      nquad = next235_cproj((2*nterms+2)*1.0d0)
c        write(*,*) nquad
      ictheta = 1
      istheta = ictheta+nquad
      icphi = istheta+nquad
      isphi = icphi+nquad
      iynm = isphi+nquad
      iynmd = iynm + (nterms+1)**2
      irat1 = iynmd + (nterms+1)**2
      irat2 = irat1 + (nterms+1)**2
      iuval = irat2 + (nterms+1)**2
      iuder = iuval + 2*nquad*(nterms+1)*(2*nterms+1)
      iephi = iuder + 2*nquad*(nterms+1)*(2*nterms+1)
      iwsave = iephi + 2*(2*nterms+1)
      iavec = iwsave + 4*nquad+20
      ibvec = iavec + 2*nquad
      lused = ibvec + 2*nquad
      if (lused.gt.lw) stop
c
      call rotmat_proj0(beta,alpha,nquad,nterms,
     1           marray2,lmpn,w(ictheta),w(istheta),
     1           w(icphi),w(isphi),w(iynm),w(iynmd),
     1           w(irat1),w(irat2),w(iuval),w(iuder),
     1           w(iephi),w(iwsave),w(iavec),w(ibvec))
      return
      end
c
c
c***********************************************************************
      subroutine rotmat_proj0(beta,alpha,nquad,nterms,
     1           marray2,lmpn,cthetas,sthetas,cphis,sphis,ynm,ynmd,
     1           rat1,rat2,uval,uder,ephis,wsave,avec,bvec)
c***********************************************************************
C
c       INPUT:
c
c       beta:  the rotation angle about the y-axis.
c       alpha:  the rotation angle about the z-axis.
c       nquad:  number of quadrature points on equator
c       nterms: order of multipole expansion
c
c       lmpn:    leading dimension of output array marray2
c       cthetas: workspace of dimension nquad
c       sthetas: workspace of dimension nquad
c       cphis:   workspace of dimension nquad
c       sphis:   workspace of dimension nquad
c       ynm:     workspace for spherical harmonics
c       ynmd:    workspace for theta derivative of spherical harmonics
c       rat1:    workspace of same dimension as ynm for precomputation
c       rat2:    workspace of same dimension as ynm for precomputation
c       uval:    workspace 
c       uder:    workspace 
c       ephis:   workspace for exp(i m phi)
c       wsave:   workspace 
c       avec:    workspace 
c       bvec:    workspace 
c
c       OUTPUT:
c
c       marray2  coefficients of rotated expansion.
c---------------------------------------------------------------------
c
      implicit real *8 (a-h,o-z)
      integer nquad,nterms
      real *8 cthetas(nquad),cphis(nquad)
      real *8 sthetas(nquad),sphis(nquad)
      real *8 ynm(0:nterms,0:nterms)
      real *8 ynmd(0:nterms,0:nterms)
      real *8 rat1(0:nterms,0:nterms)
      real *8 rat2(0:nterms,0:nterms)
      complex *16 avec(*)
      complex *16 bvec(*)
      complex *16 marray2(0:lmpn,-lmpn:lmpn)
      complex *16 uder(nquad,-nterms:nterms,0:nterms)
      complex *16 uval(nquad,-nterms:nterms,0:nterms)
      complex *16 ephis(-nterms:nterms)
      real *8 wsave(*)
c
c     Algorithm:
c     1) get locations of quadrature nodes
c     2) evaluate u and du/dtheta
c     3) project onto spherical harmonics.
c
      call getmeridian(beta,nquad,cthetas,sthetas,cphis,sphis)    
      call evalall0d(beta,alpha,nquad,cthetas,sthetas,cphis,sphis,
     2           nterms,uval,uder,ynm,ynmd,ephis,rat1,rat2)
      call projectonynm2d(nquad,
     1           uval,uder,ynm,ynmd,marray2,lmpn,nterms,
     2           m2,wsave,avec,bvec,rat1,rat2)
      return
      end
C
C
C***********************************************************************
      subroutine evalall0d(beta,alpha,nquad,cthetas,sthetas,cphis,sphis,
     1           nterms,uval,uder,ynm,ynmd,ephis,rat1,rat2)
C***********************************************************************
C
C     This subroutine evaluates the multipole expansion for each
C     order at the nquad nodes on the rotated equator.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     beta    : angle of rotation about y-axis.
c     alpha    : the rotation angle about the z-axis.
C     nquad    : number of target point son unit sphere
C     cthetas  : cos(theta) values of target points.
C     sthetas  : sin(theta) values of target points.
C     cphis    : cos(phi) values of target points.
C     sphis    : sin(phi) values of target points.
C     nterms   : order of multipole expansion
C     ynm      : work array for ynm values
C     ynmd     : work array for ynmd values
C     ephis    : work array for exp(i m phi) values
C     rat1     : work array for accelerating ynm calculation.
C     rat2     : work array for accelerating ynm calculation.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     uval(i,j) : contribution to potential 
C                 of multipole terms of order j at ith quad node.
C     uder(i,j) : contributions to theta derivative of potential
C                 of multipole terms of order j at ith quad node.
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer ndeg,morder, nquad
      real *8 cthetas(nquad),cphis(nquad)
      real *8 sthetas(nquad),sphis(nquad)
      real *8 ynm(0:nterms,0:nterms)
      real *8 ynmd(0:nterms,0:nterms)
      complex *16 ephi1,ephis(-nterms:nterms)
      complex *16 uder(nquad,-nterms:nterms,0:nterms)
      complex *16 uval(nquad,-nterms:nterms,0:nterms)
      complex *16 uv,utheta,uphi,ztmp1,ztmp2,ztsum,ztdif
      complex *16 ux,uy,uz,ima
      real *8 rat1(0:nterms,0:nterms)
      real *8 rat2(0:nterms,0:nterms)
C
      data ima/(0.0d0,1.0d0)/
      pi = 4.0d0*datan(1.0d0)
C
      cbeta = cos(beta)
      sbeta = -sin(beta)
      call ylgndrini(nterms,rat1,rat2)
c
      do i=1,nquad
      do k=0,nterms
      do j=-nterms,nterms
        uval(i,j,k)=0
        uder(i,j,k)=0
      enddo
      enddo
      enddo
c
      do jj=1,nquad/2
	 ctheta = cthetas(jj)
	 stheta = sthetas(jj)
	 cphi = cphis(jj)
	 sphi = sphis(jj)
         dir1 = -sbeta
         dir2 = 0
         dir3 = cbeta
         tang1 = cphi*ctheta
         tang2 = sphi*ctheta
         tang3 = -stheta
         proj2 = tang1*dir1 + tang2*dir2 + tang3*dir3
         tang1 = -sphi
         tang2 = cphi
         tang3 = 0
         proj1 = tang1*dir1 + tang2*dir2 + tang3*dir3
	 call ylgndru2sf(nterms,ctheta,ynm,ynmd,rat1,rat2)
         ephi1 = dcmplx(cphis(jj),sphis(jj)) *exp(-ima*alpha)
	 ephis(1) = ephi1
	 ephis(-1) = dconjg(ephi1)
	 do i = 2,nterms
	    ephis(i) = ephis(i-1)*ephi1
	    ephis(-i) = dconjg(ephis(i))
	 enddo
c
	 do ndeg = 0,nterms

	    do morder = -ndeg,ndeg

            if( morder .ne. 0 ) then
ccc            ztmp1=ephis(morder)*mpole(ndeg,morder)
            ztmp1=ephis(morder)
            uv=ynm(ndeg,abs(morder))*ztmp1*stheta
            utheta=ynmd(ndeg,abs(morder))*ztmp1
            uphi=-ynm(ndeg,abs(morder))*morder*ztmp1
            endif
            if( morder .eq. 0 ) then
ccc            uv=ynm(ndeg,0)*mpole(ndeg,0)
ccc            utheta=ynmd(ndeg,0)*stheta*mpole(ndeg,0)
            uv=ynm(ndeg,0)
            utheta=ynmd(ndeg,0)*stheta
            uphi=0
            endif
               
            irot=morder
c
c       ... apply the periodizing operator
c
            uval(jj,irot,ndeg) = uval(jj,irot,ndeg) 
     $         + uv
            uder(jj,irot,ndeg) = uder(jj,irot,ndeg) 
     $         + (utheta*proj2+uphi*ima*proj1)
            if( mod(ndeg,2) .eq. 0 ) then
            uval(jj+nquad/2,irot,ndeg) = uval(jj+nquad/2,irot,ndeg) 
     $         + uv
            uder(jj+nquad/2,irot,ndeg) = uder(jj+nquad/2,irot,ndeg) 
     $         - (utheta*proj2+uphi*ima*proj1)
            endif
            if( mod(ndeg,2) .eq. 1 ) then
            uval(jj+nquad/2,irot,ndeg) = uval(jj+nquad/2,irot,ndeg) 
     $         - uv
            uder(jj+nquad/2,irot,ndeg) = uder(jj+nquad/2,irot,ndeg)
     $         + (utheta*proj2+uphi*ima*proj1)
            endif
c
            enddo
c
	 enddo
      enddo
c
      return
      end
C
C
C
C
C
C
C***********************************************************************
      subroutine projectonynm2d(nquad,uval,uder,
     1           ynm,ynmd,marray,lmpn,nterms,m2,wsave,avec,bvec,
     $           rat1,rat2)
C***********************************************************************
C
C     This subroutine projects from values on equator for each multipole
C     order (uval, uder = dudthteta) 
C     onto spherical harmonics
C
C---------------------------------------------------------------------
C     INPUT:
C
C     nquad    : number of points on equator
C     uval     : F values on equator
C     uder     : dFdtheta values on equator
C     ynm      : work array for ynm values
C     ynmd     : work array for ynmd values
C     lmpn     : leading dim of marray (must exceed nterms)
C     nterms   : order of expansion
C     m2       : NOT IMPLEMENTED (for reduced number of degrees in 
C                expansion (second index)
C     wsave    : work array for FFT (dimension at least 4*nquad+20)
C     avec     : work array of length nquad for FFT (complex)
C     bvec     : work array of length nquad for FFT (complex)
C---------------------------------------------------------------------
C     OUTPUT:
C
C     marray   : rotated expansion 
C
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nquad, norder
      complex *16 ephi,ephi1,uval(nquad,-nterms:nterms,0:nterms)
      complex *16 uder(nquad,-nterms:nterms,0:nterms)
      complex *16 utheta,uphi,ztmp1,ztmp2
      complex *16 alpha,beta,ima
      complex *16 marray(0:lmpn,-lmpn:lmpn,-lmpn:lmpn)
      real *8 ynm(0:nterms,0:nterms)
      real *8 ynmd(0:nterms,0:nterms)
      real *8 wsave(*)
      complex *16 avec(*)
      complex *16 bvec(*)
      real *8 rat1(0:nterms,0:nterms)
      real *8 rat2(0:nterms,0:nterms)
C
      data ima/(0.0d0,1.0d0)/
c
cc      pi = 4.0d0*datan(1.0d0)
cc      theta = pi/2.0d0
cc      ctheta = cos(theta)
cc      stheta = sin(theta)
c
      ctheta = 0.0d0
      stheta = 1.0d0
      h = 1.0d0/nquad
      call ylgndru2sf(nterms,ctheta,ynm,ynmd,rat1,rat2)

      call zffti(nquad,wsave)
      do irot=-nterms,nterms
      do norder=0,nterms
         d=sqrt(2*norder+1.0d0)
         do ii = 1,nquad
            avec(ii) = uval(ii,irot,norder)*d+uder(ii,irot,norder)
         enddo
         call zfftf(nquad,avec,wsave)
         do m = -norder,norder
            if (m.ge.0)  alpha = avec(m+1)*h
            if (m.lt.0)  alpha = avec(nquad+m+1)*h
            marray(norder,m,irot) = alpha/
     1        (ynm(norder,abs(m))*d - (ynmd(norder,abs(m))))
         enddo
      enddo
      enddo

      return
      end
c
c
c
c
c
c***********************************************************************
      subroutine rotmat_projf90(beta,alpha,nterms,marray2,lmpn)
c***********************************************************************
c       Purpose:
c
c	Fast and stable algorithm for constructing the rotation operator about
c	the z-axis determined by angle alpha plus the y-axis determined
c	by angle beta. After rotation, the expansion pole is moved to
c	location (beta, alpha) in spherical coordinates (theta, phi).
c
c       The method is based on computing the induced potential and
c       its theta-derivative on the rotated equator
c       for each order (first index). The coefficients of  the rotated
c       expansion can then be obtained by FFT and projection.
c
c       There is some loss in speed over using recurrence relations 
c       but it is stable to all orders whereas the recurrence schemes 
c       are not.
c
C---------------------------------------------------------------------
c       INPUT:
c
c       beta:  the rotation angle about the y-axis.
c       alpha:  the rotation angle about the z-axis.
c       nterms: order of multipole expansion
C       lmpn    leading dim for marray2 (must exceed nterms)
c
C---------------------------------------------------------------------
c       OUTPUT:
c
c       marray2  coefficients of the rotation operator.
c
C---------------------------------------------------------------------
c
c
c
      implicit none
      integer nquad,ier,m1,m2,nterms,lmp,lmpn, next235_cproj
      integer ictheta,istheta,icphi,isphi,iynm,iynmd,irat1,irat2
      integer iuval,iuder,iephi,iwsave,iavec,ibvec,lused
      real *8 beta,alpha
      real *8, allocatable :: w(:)
      complex *16 marray2(0:lmpn,-lmpn:lmpn,-lmpn:lmpn)
      complex *16, allocatable :: cw(:)
c
        lmp=1
c
c      nquad = 2*nterms+2
c        write(*,*) nquad
      nquad = next235_cproj((2*nterms+2)*1.0d0)
c        write(*,*) nquad
      ictheta = 1
      istheta = ictheta+nquad
      icphi = istheta+nquad
      isphi = icphi+nquad
      iynm = isphi+nquad
      iynmd = iynm + (nterms+1)**2
      irat1 = iynmd + (nterms+1)**2
      irat2 = irat1 + (nterms+1)**2
      iwsave = irat2 + (nterms+1)**2
      lused = iwsave + 4*nquad+20
      allocate (w(lused), stat=ier)
      iuval = 1
      iuder = iuval + nquad*(nterms+1)*(2*nterms+1)
      iephi = iuder + nquad*(nterms+1)*(2*nterms+1)
      iavec = iephi + (2*nterms+1)
      ibvec = iavec + 2*nquad
      lused = ibvec + 2*nquad
      allocate (cw(lused), stat=ier)
      if (ier.ne.0) then 
         write(6,*) ' alloc failure in rotmat_projf90'
         stop
      endif
c
      call rotmat_proj0(beta,alpha,nquad,nterms,
     1           marray2,lmpn,w(ictheta),w(istheta),
     1           w(icphi),w(isphi),w(iynm),w(iynmd),
     1           w(irat1),w(irat2),cw(iuval),cw(iuder),
     1           cw(iephi),w(iwsave),cw(iavec),cw(ibvec))
      deallocate(w)
      return
      end
c
c
