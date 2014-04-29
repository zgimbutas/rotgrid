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
c    $Date: 2012-05-17 15:21:54 -0400 (Thu, 17 May 2012) $
c    $Revision: 2957 $
c
c
c     ROTATION VIA PROJECTION  (FORTRAN 77 and 90 VERSION).     
c
c     Multipole expansions for real-valued functions.
c
c     Requires FFT and Associated Legendre Function Libraries.
c
c     User-callable f77 routine is rotviaproj_real. 
c     User-callable f90 routine is rotviaproj_realf90.
c     The other routines are used internally.
c
c***********************************************************************
      subroutine rotviaproj_real
     1     (beta,nterms,m1,m2,mpole,lmp,marray2,lmpn,w,lw,lused)
c***********************************************************************
c       Purpose:
c
c	Fast and stable algorithm for applying rotation operator about
c	the y-axis determined by angle beta. 
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
c       If the rotation operator is to be used multiple times, and
c       memory is available, one can precompute and store the 
c       multipliers used in evalall (see below). This has not yet been
c       implemented.
c
C---------------------------------------------------------------------
c       INPUT:
c
c       beta:  the rotation angle about the y-axis.
c       nterms: order of multipole expansion
c
c       m1    : max m index for first expansion   
c               NOT IMPLEMENTED but integer argument must be supplied.
c       m2    : max m index for second expansion  
c               NOT IMPLEMENTED but integer argument must be supplied.
c
c               That is, parameters m1,m2 are currently ignored.
c
C       mpole   coefficients of original multiple expansion
C       lmp     leading dim for mpole (must exceed nterms)
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
      integer nquad,next235_rproj
      real *8 w(lw)
      complex *16 mpole(0:lmp,0:lmp)
      complex *16 marray2(0:lmpn,0:lmpn)
c
c      nquad = 2*nterms+2
c        write(*,*) nquad
      nquad = next235_rproj((2*nterms+2)*1.0d0)
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
      iuder = iuval + nquad*(nterms+1)
      iephi = iuder + nquad*(nterms+1)
      iwsave = iephi + 2*(nterms+1)
      iavec = iwsave + 4*nquad+20
      ibvec = iavec + 2*nquad
      lused = ibvec + 2*nquad
      if (lused.gt.lw) stop
c
      call rotviaproj0_real(beta,nquad,nterms,m1,m2,mpole,lmp,
     1           marray2,lmpn,w(ictheta),w(istheta),
     1           w(icphi),w(isphi),w(iynm),w(iynmd),
     1           w(irat1),w(irat2),w(iuval),w(iuder),
     1           w(iephi),w(iwsave),w(iavec),w(ibvec))
      return
      end
c
c
************************************************************************
      function next235_rproj(base)
      implicit none
      integer next235_rproj, numdiv
      real*8 base
c ----------------------------------------------------------------------
c     integer function next235_rproj returns a multiple of 2, 3, and 5
c
c     next235_rproj = 2^p 3^q 5^r >= base  where p>=1, q>=0, r>=0
************************************************************************
      next235_rproj = 2 * int(base/2d0+.9999d0)
      if (next235_rproj.le.0) next235_rproj = 2

100   numdiv = next235_rproj
      do while (numdiv/2*2 .eq. numdiv)
         numdiv = numdiv /2
      enddo
      do while (numdiv/3*3 .eq. numdiv)
         numdiv = numdiv /3
      enddo
      do while (numdiv/5*5 .eq. numdiv)
         numdiv = numdiv /5
      enddo
      if (numdiv .eq. 1) return
      next235_rproj = next235_rproj + 2
      goto 100
      end
c
c
c***********************************************************************
      subroutine rotviaproj0_real(beta,nquad,nterms,m1,m2,mpole,lmp,
     1           marray2,lmpn,cthetas,sthetas,cphis,sphis,ynm,ynmd,
     1           rat1,rat2,uval,uder,ephis,wsave,avec,bvec)
c***********************************************************************
C
c       INPUT:
c
c       beta:  the rotation angle about the y-axis.
c       nquad:  number of quadrature points on equator
c       nterms: order of multipole expansion
c
c       m1    : max m index for first expansion   
c               NOT IMPLEMENTED but integer argument must be supplied.
c       m2    : max m index for second expansion
c               NOT IMPLEMENTED but integer argument must be supplied.
c
c       mpole:   coefficients of original multiple expansion
c       lmp:     leading dimension of mpole
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
      complex *16 avec(nquad)
      complex *16 bvec(nquad)
      complex *16 mpole(0:lmp,0:lmp)
      complex *16 marray2(0:lmpn,0:lmpn)
      real *8 uder(nquad,0:nterms),uval(nquad,0:nterms)
      complex *16 ephis(0:nterms)
      real *8 wsave(4*nquad+20)
c
c     Algorithm:
c     1) get locations of quadrature nodes
c     2) evaluate u and du/dtheta
c     3) project onto spherical harmonics.
c
      call getmeridian_real(beta,nquad,cthetas,sthetas,cphis,sphis)    
      call evalall2_real(beta,nquad,cthetas,sthetas,cphis,sphis,mpole,
     2           lmp,nterms,uval,uder,ynm,ynmd,ephis,rat1,rat2)
      call projectonynm2_real
     1          (nquad,uval,uder,ynm,ynmd,marray2,lmpn,nterms,
     2           m2,wsave,avec,bvec,rat1,rat2)
      return
      end
C
C
C***********************************************************************
      subroutine getmeridian_real
     1     (beta,nquad,cthetas,sthetas,cphis,sphis)
C***********************************************************************
C     Purpose:
C
C           For a rotation of angle BETA about the y-axis, this
C           subroutine returns the NQUAD equispaced nodes on the 
C           rotated equator in the original coordinate system.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     beta  = angle of rotation
C     nquad = number of quadrature nodes in equator.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     cthetas = cos(theta) values in original coordinate system of 
C                    nquad equispaced nodes
C     sthetas = sin(theta) values in original coordinate system of 
C                    nquad equispaced nodes
C     cphis =  cos(phi) values in original coordinate system of 
C                    nquad equispaced nodes
C     sphis =  cos(phi) values in original coordinate system of 
C                    nquad equispaced nodes
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nquad
      real *8 cthetas(nquad)
      real *8 sthetas(nquad)
      real *8 cphis(nquad)
      real *8 sphis(nquad)
C
      pi = 4.0d0*datan(1.0d0)
C
      ca = cos(beta)
      sa = sin(beta)
      do i = 1,nquad
	 im1 = i-1
         phi = 2*pi*im1/nquad
	 theta = pi/2.0d0
         xp = cos(phi)*sin(theta)
         yp = sin(phi)*sin(theta)
         zp = cos(theta)
         x = ca*xp + sa*zp
         y = yp
         z = -sa*xp + ca*zp
         proj = sqrt(x**2+y**2)
	 if (proj.le.1.0d-16) then
	    cphis(i) = 1.0d0
	    sphis(i) = 0.0d0
	 else
	    cphis(i) = x/proj
	    sphis(i) = y/proj
	 endif
	 cthetas(i) = z
	 sthetas(i) = proj
      enddo
      return
      end
C
C***********************************************************************
      subroutine evalall0_real(beta,nquad,cthetas,sthetas,cphis,sphis,
     1           mpole,lmp,nterms,uval,uder,ynm,ynmd,ephis,rat1,rat2)
C***********************************************************************
C
C     This subroutine evaluates the multipole expansion for each
C     order at the nquad nodes on the rotated equator.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     beta    : angle of rotation about y-axis.
C     nquad    : number of target point son unit sphere
C     cthetas  : cos(theta) values of target points.
C     sthetas  : sin(theta) values of target points.
C     cphis    : cos(phi) values of target points.
C     sphis    : sin(phi) values of target points.
C     mpole    : original multipole expansion
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
      complex *16 mpole(0:lmp,0:lmp)
      complex *16 ephi1,ephis(0:nterms)
      real *8 uder(nquad,0:nterms),uval(nquad,0:nterms)
      real *8 uv,utheta,uphi
      complex *16 ztmp1,ztmp2,ztsum
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
      do jj=1,nquad
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
         ephi1 = dcmplx(cphis(jj),sphis(jj))
	 ephis(1) = ephi1
	 do i = 2,nterms
	    ephis(i) = ephis(i-1)*ephi1
	 enddo
c
	 do ndeg = 0,nterms
	    uv=0
	    utheta=0
	    uphi=0
	    do morder = 1,ndeg
               ztmp1=ephis(morder)*mpole(ndeg,morder)
	       uv=uv+ynm(ndeg,morder)*dble(ztmp1)
	       utheta=utheta+ynmd(ndeg,morder)*dble(ztmp1)
	       uphi=uphi-ynm(ndeg,morder)*morder*imag(ztmp1)
	    enddo
	    uv=2*stheta*uv+ynm(ndeg,0)*dble(mpole(ndeg,0))
	    utheta=2*utheta+ynmd(ndeg,0)*stheta*dble(mpole(ndeg,0))
            uphi=2*uphi
c
            uval(jj,ndeg) = uv
            uder(jj,ndeg) = (utheta*proj2-uphi*proj1)
	 enddo
      enddo
      return
      end
C
C
C
C***********************************************************************
      subroutine evalall1_real(beta,nquad,cthetas,sthetas,cphis,sphis,
     1           mpole,lmp,nterms,uval,uder,ynm,ynmd,ephis,rat1,rat2)
C***********************************************************************
C
C     This subroutine evaluates the multipole expansion for each
C     order at the nquad nodes on the rotated equator.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     beta    : angle of rotation about y-axis.
C     nquad    : number of target point son unit sphere
C     cthetas  : cos(theta) values of target points.
C     sthetas  : sin(theta) values of target points.
C     cphis    : cos(phi) values of target points.
C     sphis    : sin(phi) values of target points.
C     mpole    : original multipole expansion
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
      complex *16 mpole(0:lmp,0:lmp)
      complex *16 ephi1,ephis(0:nterms)
      real *8 uder(nquad,0:nterms),uval(nquad,0:nterms)
      real *8 uv,utheta,uphi
      complex *16 ztmp1,ztmp2,ztsum
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
         ephi1 = dcmplx(cphis(jj),sphis(jj))
	 ephis(1) = ephi1
	 do i = 2,nterms
	    ephis(i) = ephis(i-1)*ephi1
	 enddo
c
	 do ndeg = 0,nterms
	    uv=0
	    utheta=0
	    uphi=0
	    do morder = 1,ndeg
               ztmp1=ephis(morder)*mpole(ndeg,morder)
	       uv=uv+ynm(ndeg,morder)*dble(ztmp1)
	       utheta=utheta+ynmd(ndeg,morder)*dble(ztmp1)
	       uphi=uphi-ynm(ndeg,morder)*morder*imag(ztmp1)
	    enddo
	    uv=2*stheta*uv+ynm(ndeg,0)*dble(mpole(ndeg,0))
	    utheta=2*utheta+ynmd(ndeg,0)*stheta*dble(mpole(ndeg,0))
            uphi=2*uphi
c
c       ... apply the periodizing operator
c
            uval(jj,ndeg) = uv
            uder(jj,ndeg) = (utheta*proj2-uphi*proj1)
            if( mod(ndeg,2) .eq. 0 ) then
            uval(jj+nquad/2,ndeg) = +uval(jj,ndeg)
            uder(jj+nquad/2,ndeg) = -uder(jj,ndeg)
            endif
            if( mod(ndeg,2) .eq. 1 ) then
            uval(jj+nquad/2,ndeg) = -uval(jj,ndeg)
            uder(jj+nquad/2,ndeg) = +uder(jj,ndeg)
            endif
c
	 enddo
      enddo
      return
      end
C
C
C
C***********************************************************************
      subroutine evalall2_real(beta,nquad,cthetas,sthetas,cphis,sphis,
     1           mpole,lmp,nterms,uval,uder,ynm,ynmd,ephis,rat1,rat2)
C***********************************************************************
C
C     This subroutine evaluates the multipole expansion for each
C     order at the nquad nodes on the rotated equator.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     beta    : angle of rotation about y-axis.
C     nquad    : number of target point son unit sphere
C     cthetas  : cos(theta) values of target points.
C     sthetas  : sin(theta) values of target points.
C     cphis    : cos(phi) values of target points.
C     sphis    : sin(phi) values of target points.
C     mpole    : original multipole expansion
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
      complex *16 mpole(0:lmp,0:lmp)
      complex *16 ephi1,ephis(0:nterms)
      real *8 uder(nquad,0:nterms),uval(nquad,0:nterms)
      real *8 uv,utheta,uphi
      complex *16 ztmp1,ztmp2,ztsum
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
ccc        write(*,*) nquad, nquad/2, nquad/4
      nquad2=nquad/2
      if( mod(nquad2,2) .eq. 0 ) nquad4=nquad2/2+1
      if( mod(nquad2,2) .eq. 1 ) nquad4=nquad2/2+1
c
      do jj=1,nquad4
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
         ephi1 = dcmplx(cphis(jj),sphis(jj))
	 ephis(1) = ephi1
	 do i = 2,nterms
	    ephis(i) = ephis(i-1)*ephi1
	 enddo
c
	 do ndeg = 0,nterms
	    uv=0
	    utheta=0
	    uphi=0
	    do morder = 1,ndeg
               ztmp1=ephis(morder)*mpole(ndeg,morder)
	       uv=uv+ynm(ndeg,morder)*dble(ztmp1)
	       utheta=utheta+ynmd(ndeg,morder)*dble(ztmp1)
	       uphi=uphi-ynm(ndeg,morder)*morder*imag(ztmp1)
	    enddo
	    uv=2*stheta*uv+ynm(ndeg,0)*dble(mpole(ndeg,0))
	    utheta=2*utheta+ynmd(ndeg,0)*stheta*dble(mpole(ndeg,0))
            uphi=2*uphi
c
c       ... apply the periodizing operator
c
            uval(jj,ndeg) = uv
            uder(jj,ndeg) = (utheta*proj2-uphi*proj1)
            if( mod(ndeg,2) .eq. 0 ) then
            uval(jj+nquad/2,ndeg) = +uval(jj,ndeg)
            uder(jj+nquad/2,ndeg) = -uder(jj,ndeg)
            endif
            if( mod(ndeg,2) .eq. 1 ) then
            uval(jj+nquad/2,ndeg) = -uval(jj,ndeg)
            uder(jj+nquad/2,ndeg) = +uder(jj,ndeg)
            endif
c
	 enddo
         if_reflect=1
         if( jj .eq. 1 ) if_reflect=0
         if( jj .eq. nquad4 .and. mod(nquad2,2) .eq. 0 ) if_reflect=0
c
         if( if_reflect .eq. 1 ) then
         jjr=nquad2 - jj + 2
ccc         write(*,*) jj, jjr
         call ylgndr2pm_opt(nterms,ynm,ynmd)

	 ctheta = -cthetas(jj)
	 stheta = sthetas(jj)
	 cphi = -cphis(jj)
	 sphi = sphis(jj)
ccc         write(*,*) ctheta,cphi,sphi
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

         ephi1 = dcmplx(cphi,sphi)
	 ephis(1) = ephi1
	 do i = 2,nterms
	    ephis(i) = ephis(i-1)*ephi1
	 enddo
ccc         call prin2('ephis=*',ephis(-nterms),2*(2*nterms+1))
c
	 do ndeg = 0,nterms
	    uv=0
	    utheta=0
	    uphi=0
	    do morder = 1,ndeg
               ztmp1=ephis(morder)*mpole(ndeg,morder)
	       uv=uv+ynm(ndeg,morder)*dble(ztmp1)
	       utheta=utheta+ynmd(ndeg,morder)*dble(ztmp1)
	       uphi=uphi-ynm(ndeg,morder)*morder*imag(ztmp1)
	    enddo
	    uv=2*stheta*uv+ynm(ndeg,0)*dble(mpole(ndeg,0))
	    utheta=2*utheta+ynmd(ndeg,0)*stheta*dble(mpole(ndeg,0))
            uphi=2*uphi
c
c       ... apply the periodizing operator
c
            uval(jjr,ndeg) = uv
            uder(jjr,ndeg) = (utheta*proj2-uphi*proj1)
            if( mod(ndeg,2) .eq. 0 ) then
            uval(jjr+nquad/2,ndeg) = +uval(jjr,ndeg)
            uder(jjr+nquad/2,ndeg) = -uder(jjr,ndeg)
            endif
            if( mod(ndeg,2) .eq. 1 ) then
            uval(jjr+nquad/2,ndeg) = -uval(jjr,ndeg)
            uder(jjr+nquad/2,ndeg) = +uder(jjr,ndeg)
            endif
c
	 enddo
         endif
      enddo
      return
      end
C
C
C
C
C
C
C***********************************************************************
      subroutine projectonynm_real(nquad,uval,uder,
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
      complex *16 ephi,ephi1
      real *8 uval(nquad,0:1)
      real *8 uder(nquad,0:1)
      complex *16 utheta,uphi,ztmp1,ztmp2
      complex *16 alpha,beta,ima
      complex *16 marray(0:lmpn,0:lmpn)
      real *8 ynm(0:nterms,0:nterms)
      real *8 ynmd(0:nterms,0:nterms)
      real *8 wsave(4*nquad+20)
      real *8 avec(nquad)
      real *8 bvec(nquad)
      real *8 rat1(0:nterms,0:nterms)
      real *8 rat2(0:nterms,0:nterms)
C
      data ima/(0.0d0,1.0d0)/
C
cc      pi = 4.0d0*datan(1.0d0)
cc      theta = pi/2.0d0
cc      ctheta = cos(theta)
cc      stheta = sin(theta)
C
      ctheta = 0.0d0
      stheta = 1.0d0
      h = 1.0d0/nquad
      call ylgndru2sf(nterms,ctheta,ynm,ynmd,rat1,rat2)
      call dffti(nquad,wsave)
      do norder=0,nterms
	 do ii = 1,nquad
	    avec(ii) = uval(ii,norder)
	    bvec(ii) = uder(ii,norder)
         enddo
	 call dfftf(nquad,avec,wsave)
	 call dfftf(nquad,bvec,wsave)
         do m = 0,norder
	    if (m.eq.0)  alpha = avec(m+1)*h
	    if (m.eq.0)  beta = bvec(m+1)*h
	    if (m.gt.0)  alpha = dcmplx(avec(2*m),avec(2*m+1))*h
	    if (m.gt.0)  beta = dcmplx(bvec(2*m),bvec(2*m+1))*h
            marray(norder,m) = (alpha*ynm(norder,abs(m)) -
     1        beta*ynmd(norder,abs(m)))/
     1        (ynm(norder,abs(m))**2 + (ynmd(norder,abs(m)))**2)
         enddo
      enddo
      return
      end
c
c
c
c
C***********************************************************************
      subroutine projectonynm2_real(nquad,uval,uder,
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
      complex *16 ephi,ephi1
      real *8 uval(nquad,0:1)
      real *8 uder(nquad,0:1)
      complex *16 utheta,uphi,ztmp1,ztmp2
      complex *16 alpha,beta,ima
      complex *16 marray(0:lmpn,0:lmpn)
      real *8 ynm(0:nterms,0:nterms)
      real *8 ynmd(0:nterms,0:nterms)
      real *8 wsave(4*nquad+20)
      real *8 avec(nquad)
      real *8 bvec(nquad)
      real *8 rat1(0:nterms,0:nterms)
      real *8 rat2(0:nterms,0:nterms)
C
      data ima/(0.0d0,1.0d0)/
C
cc      pi = 4.0d0*datan(1.0d0)
cc      theta = pi/2.0d0
cc      ctheta = cos(theta)
cc      stheta = sin(theta)
C
      ctheta = 0.0d0
      stheta = 1.0d0
      h = 1.0d0/nquad
      call ylgndru2sf(nterms,ctheta,ynm,ynmd,rat1,rat2)
      call dffti(nquad,wsave)
      do norder=0,nterms
         d=sqrt(2*norder+1.0d0)
	 do ii = 1,nquad
	    avec(ii) = uval(ii,norder)*d + uder(ii,norder)
ccc	    bvec(ii) = uder(ii,norder)
         enddo
	 call dfftf(nquad,avec,wsave)
ccc	 call dfftf(nquad,bvec,wsave)
         do m = 0,norder
	    if (m.eq.0)  alpha = avec(m+1)*h
ccc	    if (m.eq.0)  beta = bvec(m+1)*h
	    if (m.gt.0)  alpha = dcmplx(avec(2*m),avec(2*m+1))*h
ccc	    if (m.gt.0)  beta = dcmplx(bvec(2*m),bvec(2*m+1))*h
            marray(norder,m) = alpha/
     1        (ynm(norder,abs(m))*d - (ynmd(norder,abs(m))))
         enddo
      enddo
      return
      end
c
c
c
c
c***********************************************************************
      subroutine rotviaproj_realf90(beta,nterms,m1,m2,mpole,lmp,
     1           marray2,lmpn)
c***********************************************************************
c       Purpose:
c
c	Fast and stable algorithm for applying rotation operator about
c	the y-axis determined by angle beta.
c
c       The method is based on computing the induced potential and
c       its theta-derivative on the rotated equator
c       for each order (first index). The coefficients of  the rotated
c       expansion can then be obtained by FFT and projection.
c
c       There is some loss in speed over using recurrence relations 
c       but it is stable to all orders whereas the recurrence schemes 
c       are not.
c       If the rotation operator is to be used multiple times, and
c       memory is available, one can precompute and store the 
c       multipliers used in evalall (see below). This has not yet been
c       implemented.
c
C---------------------------------------------------------------------
c       INPUT:
c
c       beta:  the rotation angle about the y-axis.
c       nterms: order of multipole expansion
C       mpole   coefficients of original multiple expansion
C       lmp     leading dim for mpole (must exceed nterms)
C       lmpn    leading dim for marray2 (must exceed nterms)
c
C---------------------------------------------------------------------
c       OUTPUT:
c
c       marray2  coefficients of rotated expansion.
c
C---------------------------------------------------------------------
c
c
c
      implicit none
      integer nquad,ier,m1,m2,nterms,lmp,lmpn, next235_rproj
      integer ictheta,istheta,icphi,isphi,iynm,iynmd,irat1,irat2
      integer iuval,iuder,iephi,iwsave,iavec,ibvec,lused
      real *8 beta
      real *8, allocatable :: w(:)
      complex *16 mpole(0:lmp,0:lmp)
      complex *16 marray2(0:lmpn,0:lmpn)
      complex *16, allocatable :: cw(:)
c
c      nquad = 2*nterms+2
c        write(*,*) nquad
      nquad = next235_rproj((2*nterms+2)*1.0d0)
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
      iuder = iuval + nquad*(nterms+1)
      iephi = iuder + nquad*(nterms+1)
      iavec = iephi + (2*nterms+1)
      ibvec = iavec + 2*nquad
      lused = ibvec + 2*nquad
      allocate (cw(lused), stat=ier)
      if (ier.ne.0) then 
         write(6,*) ' alloc failure in rotviaprojf90'
         stop
      endif
c
      call rotviaproj0_real
     $          (beta,nquad,nterms,nterms,nterms,mpole,lmp,
     1           marray2,lmpn,w(ictheta),w(istheta),
     1           w(icphi),w(isphi),w(iynm),w(iynmd),
     1           w(irat1),w(irat2),cw(iuval),cw(iuder),
     1           cw(iephi),w(iwsave),cw(iavec),cw(ibvec))
      return
      end
c
c
