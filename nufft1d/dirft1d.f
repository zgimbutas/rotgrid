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
c***********************************************************************
      subroutine dirft1d1(nj,xj,cj, iflag, ms,fk)
      implicit none
      integer nj, iflag, ms
      real*8 xj(nj)
      complex*16 cj(nj), fk(-ms/2:(ms-1)/2)
c ----------------------------------------------------------------------
c     direct computation of nonuniform FFT
c
c                  1  nj
c     fk(k1)    = -- SUM cj(j) exp(+/-i k1 xj(j)) 
c                 nj j=1
c
c     for -ms/2 <= k1 <= (ms-1)/2
c
c     If (iflag .ge.0) the + sign is used in the exponential.
c     If (iflag .lt.0) the - sign is used in the exponential.
c
c***********************************************************************
      integer j, k1
      complex*16 zf, cm1
c
      do k1 = -ms/2, (ms-1)/2
         fk(k1) = dcmplx(0d0,0d0)
      enddo
c
      do j = 1, nj
         if (iflag .ge. 0) then
            zf = dcmplx(dcos(xj(j)),+dsin(xj(j)))
         else
            zf = dcmplx(dcos(xj(j)),-dsin(xj(j)))
         endif
c
         cm1 = cj(j) / dble(nj)
         do k1 = 0, (ms-1)/2
            fk(k1) = fk(k1) + cm1
            cm1 = cm1 * zf
         enddo
c
         zf = dconjg(zf)
         cm1 = cj(j) / dble(nj)
         do k1 = -1, -ms/2, -1
            cm1 = cm1 * zf
            fk(k1) = fk(k1) + cm1
         enddo
      enddo
      end
c
c
c
c
c
c***********************************************************************
      subroutine dirft1d2(nj,xj,cj, iflag, ms,fk)
      implicit none
      integer nj, iflag, ms
      real*8 xj(nj)
      complex*16 cj(nj), fk(-ms/2:(ms-1)/2)
c ----------------------------------------------------------------------
c     direct computation of nonuniform FFT
c
c     cj(j) = SUM   fk(k1) exp(+/-i k1 xj(j)) 
c             k1  
c                            for j = 1,...,nj
c
c     where -ms/2 <= k1 <= (ms-1)/2
c
c
c     If (iflag .ge.0) the + sign is used in the exponential.
c     If (iflag .lt.0) the - sign is used in the exponential.
c***********************************************************************
      integer j, k1
      complex*16 zf, cm1
c
      do j = 1, nj
         if (iflag .ge. 0) then
            zf = dcmplx(dcos(xj(j)),+dsin(xj(j)))
         else
            zf = dcmplx(dcos(xj(j)),-dsin(xj(j)))
         endif
c
         cj(j) = fk(0)
         cm1 = zf
         do k1 = 1, (ms-1)/2
            cj(j) = cj(j) + cm1*fk(k1)+dconjg(cm1)*fk(-k1)
            cm1 = cm1 * zf
         enddo
         if (ms/2*2.eq.ms) cj(j) = cj(j) + dconjg(cm1)*fk(-ms/2)
      enddo
      end
c
c
c
c
c
c
c***********************************************************************
      subroutine dirft1d3(nj,xj,cj, iflag, nk,sk,fk)
      implicit none
      integer nj, iflag, nk
      real*8 xj(nj), sk(nk)
      complex*16 cj(nj), fk(nk)
c ----------------------------------------------------------------------
c     direct computation of nonuniform FFT
c
c              nj
c     fk(k) = SUM cj(j) exp(+/-i s(k) xj(j)) 
c             j=1                   
c
c                    for k = 1, ..., nk
c
c     If (iflag .ge.0) the + sign is used in the exponential.
c     If (iflag .lt.0) the - sign is used in the exponential.
c
c***********************************************************************
      integer j, k
      real*8 ssk
c
      do k = 1, nk
         if (iflag .ge. 0) then
            ssk =  sk(k)
         else
            ssk =  -sk(k)
         endif
c
         fk(k) = dcmplx(0d0, 0d0)
         do j = 1, nj
            fk(k) = fk(k) + cj(j)*dcmplx(cos(ssk*xj(j)),sin(ssk*xj(j)))
         enddo
      enddo
      end
