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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c
c
c     ROTATION VIA RECURSION  (FORTRAN 77 AND 90 VERSIONS).
c     good up to order 100 or so.
c
c     Multipole expansions for real-valued functions.
c
c       rotviarecur3pa_apply_real - apply the rotation matrix, move pole to 
c                              location (beta, alpha)
c       rotviarecur3pf_apply_real - apply the rotation matrix, move pole to 
c                        locations (beta, alpha_k), alpha_k = 2*pi *(k-1)/nrot
c
C*****************************************************************
        subroutine rotviarecur3pa_apply_real
     $     (theta,alpha,nterms,m1,m2,mpole,ld1,marray,ld2,rotmat,ldc)
C*****************************************************************
c
c       Purpose:
c
c	Fast, recursive method for applying rotation matrix about
c       the z-axis determined by angle alpha plus the y-axis determined
c       by angle theta. After rotation, the expansion pole is moved to
c       location (theta, alpha) in spherical coordinates.
c
C---------------------------------------------------------------------
c       INPUT:
c
c       theta:  the rotation angle about the y-axis.
c       alpha:  the rotation angle about the z-axis.
c       nterms: order of multipole expansion
c
c       m1    :  max m index for first expansion  
c       m2    :  max m index for second expansion 
C       mpole :  coefficients of original multipole expansion
C       ld1   :  leading dim for mpole (must exceed nterms)
C       ld2   :  leading dim for marray (must exceed nterms)
c       rotmat:  real *8 (0:ldc,0:ldc,-ldc:ldc): rotation matrix 
c       ldc   :  leading dim for rotation matrix (must exceed nterms)
c
C---------------------------------------------------------------------
c       OUTPUT:
c
c       marray   coefficients of rotated expansion.
c
C---------------------------------------------------------------------
	implicit none
	integer nterms,m1,m2,ld1,ld2,ldc
	integer ij,m,mp
	real *8 theta,alpha
	complex *16 mpole(0:ld1,0:ld1)
	complex *16 marray(0:ld2,0:ld2)
	real *8 rotmat(0:ldc,0:ldc,-ldc:ldc)
        complex *16 ima,cd,cd0
        data ima/(0.0d0,1.0d0)/
c
        if (m1.lt.nterms .or. m2.lt.nterms) then
c
c       ... truncate multipole expansions
c
           do ij=0,nterms

              do m=0,ij
                  marray(ij,m)=0
              enddo
c
              do m=0,min(ij,m2)
                 marray(ij,m)=mpole(ij,0)*rotmat(0,ij,m) 
                 do mp=1,min(ij,m1)
                    cd=exp(+ima*mp*alpha)
                    marray(ij,m)=marray(ij,m)+
     1	            mpole(ij,mp)*cd*rotmat(mp,ij,m)+
     1              dconjg(mpole(ij,mp))*dconjg(cd)*rotmat(mp,ij,-m)
                 enddo
              enddo
           enddo
        else
c
c       ... apply rotation matrix
c
           do ij=0,nterms
c
c              do m=0,ij
c                 marray(ij,m)=0
c              enddo
c
              do m=0,ij
                 marray(ij,m)=mpole(ij,0)*rotmat(0,ij,m) 
                 cd0=exp(+ima*alpha)
                 cd=cd0
                 do mp=1,ij
                    marray(ij,m)=marray(ij,m)+
     1	            mpole(ij,mp)*cd*rotmat(mp,ij,m)+
     1              dconjg(mpole(ij,mp))*dconjg(cd)*rotmat(mp,ij,-m)
                    cd=cd*cd0
                 enddo
              enddo
           enddo
        endif
        return
        end
c        
c
c
c
c
C*****************************************************************
        subroutine rotviarecur3px_apply_real
     $     (theta,nterms,m1,m2,mpole,ld1,marray,ld2,rotmat,ldc)
C*****************************************************************
c
c       Purpose:
c
c	Fast, recursive method for applying rotation matrix about
c	the y-axis determined by angle theta.
c
C---------------------------------------------------------------------
c       INPUT:
c
c       theta:  the rotation angle about the y-axis.
c       nterms: order of multipole expansion
c
c       m1    :  max m index for first expansion  
c       m2    :  max m index for second expansion 
C       mpole :  coefficients of original multipole expansion
C       ld1   :  leading dim for mpole (must exceed nterms)
C       ld2   :  leading dim for marray (must exceed nterms)
c       rotmat:  real *8 (0:ldc,0:ldc,-ldc:ldc): rotation matrix 
c       ldc   :  leading dim for rotation matrix (must exceed nterms)
c
C---------------------------------------------------------------------
c       OUTPUT:
c
c       marray   coefficients of rotated expansion.
c
C---------------------------------------------------------------------
	implicit none
	integer nterms,m1,m2,ld1,ld2,ldc
	integer ij,m,mp
	real *8 theta
	complex *16 mpole(0:ld1,0:ld1)
	complex *16 marray(0:ld2,0:ld2)
	real *8 rotmat(0:ldc,-ldc:ldc,-ldc:ldc)
c
        if (m1.lt.nterms .or. m2.lt.nterms) then
c
c       ... truncate multipole expansions
c
           do ij=0,nterms

              do m=0,ij
                  marray(ij,m)=0
              enddo
c
              do m=0,min(ij,m2)
                 marray(ij,m)=mpole(ij,0)*rotmat(ij,m,0) 
                 do mp=1,min(ij,m1)
                    marray(ij,m)=marray(ij,m)+
     1	            mpole(ij,mp)*rotmat(ij,m,mp)+
     1              dconjg(mpole(ij,mp))*rotmat(ij,-m,mp)
                 enddo
              enddo
           enddo
        else
c
c       ... apply rotation matrix
c
           do ij=0,nterms
c
c              do m=0,ij
c                 marray(ij,m)=0
c              enddo
c
c              do m=0,ij
c                 marray(ij,m)=mpole(ij,0)*rotmat(ij,m,0) 
c                 do mp=1,ij
c                    marray(ij,m)=marray(ij,m)+
c     1	             mpole(ij,mp)*rotmat(ij,m,mp)+
c     1              dconjg(mpole(ij,mp))*rotmat(ij,-m,mp)
ccc                    write(*,*) rotmat(ij,-m,mp), rotmat(ij,m,-mp)
c                 enddo
c              enddo
c
              do m=0,ij
                 marray(ij,m)=mpole(ij,0)*rotmat(ij,m,0) 
              enddo
              do mp=1,ij
                 do m=0,ij
                    marray(ij,m)=marray(ij,m)+
     1	            mpole(ij,mp)*rotmat(ij,m,mp)+
     1              dconjg(mpole(ij,mp))*rotmat(ij,-m,mp)
                 enddo
              enddo
c
c
           enddo
        endif
        return
        end
c        
c
c
C*****************************************************************
        subroutine rotviarecur3pf_apply_real
     $     (theta,nrot,nterms,
     $     m1,m2,mpole,ld1,marray,ld2,rotmat,ldc)
C*****************************************************************
c
c       Purpose:
c
c       Fast, recursive algorithm for applying rotation operators about
c       the z-axis determined by angles alpha_k plus the y-axis determined
c       by angle beta. After rotation, the expansion pole is moved to nrot
c       locations (beta, alpha_k) in spherical coordinates (theta, phi),
c       where alpha_k = 2*pi * k/nrot, k=0..nrot-1.
c
C---------------------------------------------------------------------
c       INPUT:
c
c       theta:  the rotation angle about the y-axis.
c       nrot:  the number of the rotation angles about the z-axis.
c       nterms: order of multipole expansion
c
c       m1    :  max m index for first expansion  
c       m2    :  max m index for second expansion 
C       mpole :  coefficients of original multipole expansion
C       ld1   :  leading dim for mpole (must exceed nterms)
C       ld2   :  leading dim for marray (must exceed nterms)
c       rotmat:  real *8 (0:ldc,0:ldc,-ldc:ldc): rotation matrix 
c       ldc   :  leading dim for rotation matrix (must exceed nterms)
c
C---------------------------------------------------------------------
c       OUTPUT:
c
c       marray   coefficients of rotated expansions.
c
C---------------------------------------------------------------------
ccc	implicit none
        implicit real *8 (a-h,o-z)
	integer nterms,m1,m2,ld1,ld2,ldc
	integer ij,m,mp,nrot
	real *8 theta,alpha,pi,done
	complex *16 mpole(0:ld1,0:ld1)
	complex *16 marray(0:ld2,0:ld2,nrot)
	real *8 rotmat(0:ldc,0:ldc,-ldc:ldc)
        complex *16 ima,ztmp,cd
        real *8, allocatable :: wsave(:)
        complex *16, allocatable :: w(:)
        data ima/(0.0d0,1.0d0)/
c
        done=1
        pi=4*atan(done)

        allocate( wsave(4*nrot+15) )
        call zffti(nrot,wsave)

        allocate( w(nrot) )
c

        if( nrot .lt. 4 ) then
c
c       ... apply rotation matrices directly
c
        do k=1,nrot
c
        alpha = 2*pi*(k-1)/dble(nrot)
c
        do ij=0,nterms
c
        do m=0,ij
          marray(ij,m,k)=mpole(ij,0)*rotmat(0,ij,m) 
          do mp=1,ij
          cd=exp(+ima*mp*alpha)
          marray(ij,m,k)=marray(ij,m,k)+
     1       mpole(ij,mp)*cd*rotmat(mp,ij,m)+
     1       dconjg(mpole(ij,mp))*dconjg(cd)*rotmat(mp,ij,-m)
          enddo
        enddo

        enddo

        enddo

        endif
c
c       
        if( nrot .ge. 4 .and. nrot .lt. 2*nterms+1 ) then
c
c       ... apply rotation matrices via FFT 
c       with periodization (aliasing)
c
        do ij=0,nterms
c
        do m=0,ij

        do mp = 1,nrot
          w(mp)=0
        enddo
c
        w(1) = mpole(ij,0)*rotmat(0,ij,m)
        mp1 = 2
        mp2 = nrot
        do mp=1,ij          
        w(mp1) = w(mp1)+mpole(ij,mp)*rotmat(mp,ij,m) 
        w(mp2) = w(mp2)+dconjg(mpole(ij,mp))*rotmat(mp,ij,-m) 
        mp1=mp1+1
        if( mp1 .gt. nrot ) mp1=1
        mp2=mp2-1
        if( mp2 .lt. 1 ) mp2=nrot
        enddo
        call zfftb(nrot,w,wsave)

        do k=1,nrot
        marray(ij,m,k) = w(k)
        enddo

        enddo
        
        enddo
c
        endif
c
c
c
        if( nrot .ge. 2*nterms+1 ) then
c
c       ... apply rotation matrices via FFT
c       with zero padding
c
        do ij=0,nterms
c
        do m=0,ij

        do mp = 1,nrot
          w(mp)=0
        enddo
c
        w(1) = mpole(ij,0)*rotmat(0,ij,m)
        do mp=1,ij          
          w(mp+1) = mpole(ij,mp)*rotmat(mp,ij,m) 
          w(nrot-mp+1) = dconjg(mpole(ij,mp))*rotmat(mp,ij,-m) 
        enddo
        call zfftb(nrot,w,wsave)

        do k=1,nrot
        marray(ij,m,k) = w(k)
        enddo

        enddo
        
        enddo
c
        endif
c
        return
        end
c        
c
c
c
c
C*****************************************************************
        subroutine rot1lat_wfft_real
     $     (theta,nrot,nterms,
     $     m1,m2,mpole,ld1,marray,ld2,rotmat,ldc)
C*****************************************************************
c
c       Purpose:
c
c       Fast, recursive algorithm for applying rotation operators about
c       the z-axis determined by angles alpha_k plus the y-axis determined
c       by angle beta. After rotation, the expansion pole is moved to nrot
c       locations (beta, alpha_k) in spherical coordinates (theta, phi),
c       where alpha_k = 2*pi * k/nrot, k=0..nrot-1.
c
c       NOTE: this routine uses different memory layout for rotmat,
c       generated by rotmat_proj routine.
c
C---------------------------------------------------------------------
c       INPUT:
c
c       theta:  the rotation angle about the y-axis.
c       nrot:  the number of the rotation angles about the z-axis.
c       nterms: order of multipole expansion
c
c       m1    :  max m index for first expansion  
c       m2    :  max m index for second expansion 
C       mpole :  coefficients of original multipole expansion
C       ld1   :  leading dim for mpole (must exceed nterms)
C       ld2   :  leading dim for marray (must exceed nterms)
c       rotmat:  real *8 (0:ldc,-ldc:ldc,-ldc:ldc): rotation matrix 
c       ldc   :  leading dim for rotation matrix (must exceed nterms)
c
C---------------------------------------------------------------------
c       OUTPUT:
c
c       marray   coefficients of rotated expansions.
c
C---------------------------------------------------------------------
ccc	implicit none
        implicit real *8 (a-h,o-z)
	integer nterms,m1,m2,ld1,ld2,ldc
	integer ij,m,mp,nrot
	real *8 theta,alpha,pi,done
	complex *16 mpole(0:ld1,0:ld1)
	complex *16 marray(0:ld2,0:ld2,nrot)
	real *8 rotmat(0:ldc,-ldc:ldc,-ldc:ldc)
        complex *16 ima,ztmp,cd
        real *8, allocatable :: wsave(:)
        complex *16, allocatable :: w(:)
        data ima/(0.0d0,1.0d0)/
c
        done=1
        pi=4*atan(done)

        allocate( wsave(4*nrot+15) )
        call zffti(nrot,wsave)

        allocate( w(nrot) )
c

        if( nrot .lt. 4 ) then
c
c       ... apply rotation matrices directly
c
        do k=1,nrot
c
        alpha = 2*pi*(k-1)/dble(nrot)
c
        do ij=0,nterms
c
        do m=0,ij
          marray(ij,m,k)=mpole(ij,0)*rotmat(ij,m,0) 
          do mp=1,ij
          cd=exp(+ima*mp*alpha)
          marray(ij,m,k)=marray(ij,m,k)+
     1       mpole(ij,mp)*cd*rotmat(ij,m,mp)+
     1       dconjg(mpole(ij,mp))*dconjg(cd)*rotmat(ij,-m,mp)
          enddo
        enddo

        enddo

        enddo

        endif
c
c       
        if( nrot .ge. 4 .and. nrot .lt. 2*nterms+1 ) then
c
c       ... apply rotation matrices via FFT 
c       with periodization (aliasing)
c
        do ij=0,nterms
c
        do m=0,ij

        do mp = 1,nrot
          w(mp)=0
        enddo
c
        w(1) = mpole(ij,0)*rotmat(ij,m,0)
        mp1 = 2
        mp2 = nrot
        do mp=1,ij          
        w(mp1) = w(mp1)+mpole(ij,mp)*rotmat(ij,m,mp) 
        w(mp2) = w(mp2)+dconjg(mpole(ij,mp))*rotmat(ij,-m,mp) 
        mp1=mp1+1
        if( mp1 .gt. nrot ) mp1=1
        mp2=mp2-1
        if( mp2 .lt. 1 ) mp2=nrot
        enddo
        call zfftb(nrot,w,wsave)

        do k=1,nrot
        marray(ij,m,k) = w(k)
        enddo

        enddo
        
        enddo
c
        endif
c
c
c
        if( nrot .ge. 2*nterms+1 ) then
c
c       ... apply rotation matrices via FFT
c       with zero padding
c
        do ij=0,nterms
c
        do m=0,ij

        do mp = 1,nrot
          w(mp)=0
        enddo
c
        w(1) = mpole(ij,0)*rotmat(ij,m,0)
        do mp=1,ij          
          w(mp+1) = mpole(ij,mp)*rotmat(ij,m,mp) 
          w(nrot-mp+1) = dconjg(mpole(ij,mp))*rotmat(ij,-m,mp) 
        enddo
        call zfftb(nrot,w,wsave)

        do k=1,nrot
        marray(ij,m,k) = w(k)
        enddo

        enddo
        
        enddo
c
        endif
c
        return
        end
c        
c
c
c
c
