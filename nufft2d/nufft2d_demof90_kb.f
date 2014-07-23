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
      program testfft
ccc      implicit none
      implicit real *8 (a-h,o-z)
c
      integer i,ier,iflag,j,k1,k2,mx,ms,mt,n1,n2,nj,nk
ccc      parameter (mx=2560*2560/4)
ccc      parameter (mx=2560*2560)
      parameter (mx=5476*6400)
      real*8 xj(mx),yj(mx)
      real *8 sk(mx),tk(mx)
      real*8 err,pi,eps,salg,ealg
      parameter (pi=3.141592653589793238462643383279502884197d0)

c      complex*16 cj(mx),cj0(mx),cj1(mx)
c      complex*16 fk0(mx),fk1(mx)

c
c      real *8 fm1(mx), xc(21*mx)
c      real *8 fm2(mx), yc(21*mx)

      complex *16, allocatable :: cj(:), cj0(:), cj1(:)
      complex *16, allocatable :: fk0(:), fk1(:)

      real *8, allocatable :: w(:)

        allocate(cj(mx))
        allocate(cj0(mx))
        allocate(cj1(mx))

        allocate(fk0(mx))
        allocate(fk1(mx))

        allocate(w(2*mx+2*15*mx))
        lw = 2*mx+2*15*mx

c     --------------------------------------------------
c     create some test data
c     --------------------------------------------------
c
      n1 = 360*2
      n2 = 400*2
      ms = 32
      mt = 30

c      n1 = 30**2
c      n2 = 30**2
c      ms = 30
c      mt = 30


      np=24+1
      ms = np*2
      mt = np*2+1
      mt = next235(mt*1.0d0)
      n1 = ms**2/4
      n2 = mt**2

      nj = n1*n2
      write(*,*) ms,mt,nj
c
      do k1 = -n1/2, (n1-1)/2
         do k2 = -n2/2, (n2-1)/2
            j = (k2+n2/2+1) + (k1+n1/2)*n2
            xj(j) = pi*dcos(-pi*k1/n1)
            yj(j) = pi*dcos(-pi*k2/n2)
            cj(j) = dcmplx(dsin(pi*j/n1),dcos(pi*j/n2))
         enddo
      enddo
c
c     -----------------------
c     start tests
c     -----------------------
c
      iflag = 1
      print*,'Starting 2D testing: ', ' nj =',nj, ' ms,mt =',ms,mt

         call dirft2d1(nj,xj,yj,cj,iflag,ms,mt,fk0)
         call dirft2d2(nj,xj,yj,cj0,iflag,ms,mt,fk0)

      do i = 1,4
         if (i.eq.1) eps=1d-4
         if (i.eq.2) eps=1d-8
         if (i.eq.3) eps=1d-12
         if (i.eq.4) eps=1d-16
c extented/quad precision tests
         if (i.eq.5) eps=1d-20
         if (i.eq.6) eps=1d-24
         if (i.eq.7) eps=1d-28
         if (i.eq.8) eps=1d-32
	 print*,' '
	 print*,' Requested precision eps =',eps
	 print*,' '
c
c     -----------------------
c     call 2D Type 1 method
c     -----------------------
c

c
         t1=second()
         call nufft2dkbi(nj,xj,yj,eps,ms,mt,ier,w,lw,lused)
         t2=second()
         write(*,*) 'init=*',t2-t1
         
         t1=second()
         call nufft2d1f90(nj,xj,yj,cj,iflag,eps,ms,mt,fk1,ier)
         t2=second()
         write(*,*) 'time=*',t2-t1

         call errcomp(fk0,fk1,ms*mt,err)
         print *, ' ier = ',ier
         call errcomp(fk0,fk1,ms*mt,err)
         print *, ' type 1 err = ',err

         t1=second()
         call nufft2d1f90kbw(nj,xj,yj,cj,iflag,eps,ms,mt,fk1,ier,w)
         t2=second()
         write(*,*) 'time=*',t2-t1

         call errcomp(fk0,fk1,ms*mt,err)
         print *, ' ier = ',ier
         call errcomp(fk0,fk1,ms*mt,err)
         print *, ' type 1 err = ',err

c
c
c     -----------------------
c      call 2D Type 2 method
c     -----------------------

         t1=second()
         call nufft2d2f90(nj,xj,yj,cj1,iflag,eps,ms,mt,fk0,ier)
         t2=second()
         write(*,*) 'time=*',t2-t1

         print *, ' ier = ',ier
         call errcomp(cj0,cj1,nj,err)
         print *, ' type 2 err = ',err

         t1=second()
         call nufft2d2f90kbw(nj,xj,yj,cj1,iflag,eps,ms,mt,fk0,ier,w)
         t2=second()
         write(*,*) 'time=*',t2-t1

         print *, ' ier = ',ier
         call errcomp(cj0,cj1,nj,err)
         print *, ' type 2 err = ',err

         cycle
c
c     -----------------------
c      call 2D Type3 method
c     -----------------------
         nk = ms*mt
         do k1 = 1, nk
            sk(k1) = 48*(dcos(k1*pi/nk))
            tk(k1) = 32*(dsin(-pi/2+k1*pi/nk))
         enddo

         call dirft2d3(nj,xj,yj,cj,iflag,nk,sk,tk,fk0)
         call nufft2d3f90(nj,xj,yj,cj,iflag,eps,nk,sk,tk,fk1,ier)
c
         print *, ' ier = ',ier
         call errcomp(fk0,fk1,nk,err)
         print *, ' type 1 err = ',err
      enddo 
      stop
      end
c
c
c
c
c
      subroutine errcomp(fk0,fk1,n,err)
      implicit none
      integer k,n
      complex*16 fk0(n), fk1(n)
      real *8 salg,ealg,err
c
      ealg = 0d0
      salg = 0d0
      do k = 1, n
         ealg = ealg + cdabs(fk1(k)-fk0(k))**2
         salg = salg + cdabs(fk0(k))**2
      enddo
      err =sqrt(ealg/salg)
      return
      end
