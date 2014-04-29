cc Copyright (C) 2011-2012: Zydrunas Gimbutas
cc Contact: gimbutas@cims.nyu.edu
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
c
c
        subroutine rotgrid_fsr_cmpl_init(nterms,nbeta,beta,rotmat)
        implicit real *8 (a-h,o-z)
c
c    Initialize rotgrid_fsr_cmpl.
c  
c    Fast, FFT-based algorithm for rotating complex spherical harmonic grids
c    into pole locations (beta_j, alpha_k), 
c    where alpha_k = 2*pi * k/nrot, k=0..nrot-1, beta_j=1..nbeta
c  
c    Input parameters:
c  
c    nterms - the number of terms in spherical harmonics expansion
c    nbeta - the number of points in pole location meridian discretization
c    beta - angles for new pole locations beta_j, j=1..nbeta
c  
c    Output parameters:
c  
c    rotmat - The rotation operators for directions (beta_j,0). 
c
        real *8 beta(nbeta)
        real *8 rotmat(0:nterms,-nterms:nterms,-nterms:nterms,nbeta) 
        complex *16, allocatable :: rotmat_cmpl(:,:,:)
        
        allocate( rotmat_cmpl(0:nterms,-nterms:nterms,-nterms:nterms) )
c
        alpha=0
c
        do kk=1,nbeta
        call rotmat_projf90(beta(kk),alpha,nterms,rotmat_cmpl,nterms)
c
        do j=-nterms,nterms
        do i=-nterms,nterms
        do k=0,nterms
        rotmat(k,i,j,kk) = rotmat_cmpl(k,i,j)
        enddo
        enddo
        enddo

        enddo
c
        return
        end
c
c
c
c
c
        subroutine rotgrid_fsr_cmpl(nterms,mpole,nphi,ntheta,
     $     nbeta,beta,nrot,grids,rotmat,ctheta,ynms,wsave)
        implicit real *8 (a-h,o-z)
c
c   Rotate the complex spherical harmonic grids.
c
c  Fast, FFT-based algorithm for rotating complex spherical harmonic grids
c  (dimensioned NPHI-by-NTHETA) into new pole locations (beta_j, alpha_k), 
c  where alpha_k = 2*pi * k/nrot, k=0..nrot-1,  beta_j, j=1..nbeta. 
c
c  GRIDS = rotgrid_fsr_cmpl(NTERMS,MPOLE,NPHI,NTHETA,NBETA,BETA,NROT,...
c       ROTMAT,CTHETA,YNMS,WSAVE) rotates 
c  the complex spherical harmonics expansion of degree NTERMS 
c  about the z-axis by degree ALPHA_K and about the y-axis by degree BETA.
c  into a collection of new pole locations (beta_j, alpha_k) in spherical 
c  coordinates (theta, phi), where alpha_k = 2*pi * k/nrot, k=0..nrot-1,
c  beta_j, j=1..nbeta.
c
c  The rotated poles form a uniformly spaced grid on lattitude \theta.
c
c  grids - function values on the rotated grids, 
c             NPHI-by-NTHETA-by-NROT-by-NBETA complex*16 matrix
c
c       Input parameters:
c
c  nterms - the number of terms in spherical harmonics expansion
c  mpole - the coefficients of spherical harmonics expansion,
c                    complex*16 (0:nterms,-nterms:nterms)
c  nphi - the number of points in latitude discretization (for spherical grid)
c  ntheta - the number of points in meridian discretization (for spherical grid)
c  ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
c  ynms - Legendre functions, real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
c  wsave - fftpack wsave array, complex *16 wsave(4*nphi+15)
c
c  nbeta - the number of points in pole location meridian discretization
c  beta - angles for new pole locations beta_j, j=1..nbeta
c  nrot - angles for new pole locations alpha_k = 2*pi * k/nrot, k=0..nrot-1.
c       
c  rotmat - The rotation operators for directions (beta_j,0). 
c
c  [rotmat, ctheta, ynms, wsave] must be initialized via 
c              a preceding call to rotgrid_fsr_cmpl_init
c
c       Output parameters:
c
c  grids - function values on the rotated grids, 
c             NPHI-by-NTHETA-by-NROT-by-NBETA complex*16 matrix
c
c  Our definition of complex spherical harmonics is
c
c  Ynm(theta,phi)= sqrt( 2n+1) sqrt((n-m)!/(n+m)!) 
c                  Pnm(cos theta) e^(im phi), 
c  Yn,-m(theta,phi) = sqrt( 2n+1) sqrt((n-m)!/(n+m)!) 
c                  Pnm(cos theta) e^(-im phi),   for m >= 0.
c       
c  Note that we do not include the Condon-Shortley phase (-1)^m, if m<0.
c
        complex *16 mpole(0:nterms,-nterms:nterms)

        real *8 beta(nbeta)
        real *8 rotmat(0:nterms,-nterms:nterms,-nterms:nterms,nbeta) 
        real *8 ctheta(ntheta),ynms(0:nterms,0:nterms,ntheta/2+1)
        complex *16 wsave(4*nphi+15)
c
        complex *16 grids(nphi,ntheta,nrot,nbeta)
c
        complex *16, allocatable :: marray(:,:,:) 
c
        allocate( marray(0:nterms,-nterms:nterms,nrot) )
c
        do k=1,nbeta
        call rot1lat_wfft
     $     (beta(k),nrot,nterms,
     $     nterms,nterms,mpole,nterms,marray,
     $     nterms,rotmat(0,-nterms,-nterms,k),nterms)
        do i=1,nrot
        call sphtrans_cmpl(
     $     nterms,marray(0,-nterms,i),
     $     nphi,ntheta,grids(1,1,i,k),
     $     ctheta,ynms,wsave)
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine rotgrid_dsr_cmpl(nterms,mpole,nphi,ntheta,
     $     nbeta,beta,nrot,grids,rotmat,ctheta,ynms,wsave)
        implicit real *8 (a-h,o-z)
c
c   Rotate the complex spherical harmonic grids.
c
c  Fast, direct algorithm for rotating complex spherical harmonic grids
c  (dimensioned NPHI-by-NTHETA) into new pole locations (beta_j, alpha_k), 
c  where alpha_k = 2*pi * k/nrot, k=0..nrot-1,  beta_j, j=1..nbeta. 
c
c  GRIDS = rotgrid_fsr_cmpl(NTERMS,MPOLE,NPHI,NTHETA,NBETA,BETA,NROT,...
c       ROTMAT,CTHETA,YNMS,WSAVE) rotates 
c  the complex spherical harmonics expansion of degree NTERMS 
c  about the z-axis by degree ALPHA_K and about the y-axis by degree BETA.
c  into a collection of new pole locations (beta_j, alpha_k) in spherical 
c  coordinates (theta, phi), where alpha_k = 2*pi * k/nrot, k=0..nrot-1,
c  beta_j, j=1..nbeta.
c
c  The rotated poles form a uniformly spaced grid on lattitude \theta.
c
c  grids - function values on the rotated grids, 
c             NPHI-by-NTHETA-by-NROT-by-NBETA complex*16 matrix
c
c       Input parameters:
c
c  nterms - the number of terms in spherical harmonics expansion
c  mpole - the coefficients of spherical harmonics expansion,
c                    complex*16 (0:nterms,-nterms:nterms)
c  nphi - the number of points in latitude discretization (for spherical grid)
c  ntheta - the number of points in meridian discretization (for spherical grid)
c  ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
c  ynms - Legendre functions, real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
c  wsave - fftpack wsave array, complex *16 wsave(4*nphi+15)
c
c  nbeta - the number of points in pole location meridian discretization
c  beta - angles for new pole locations beta_j, j=1..nbeta
c  nrot - angles for new pole locations alpha_k = 2*pi * k/nrot, k=0..nrot-1.
c       
c  rotmat - The rotation operators for directions (beta_j,0). 
c
c  [rotmat, ctheta, ynms, wsave] must be initialized via 
c              a preceding call to rotgrid_fsr_cmpl_init
c
c       Output parameters:
c
c  grids - function values on the rotated grids, 
c             NPHI-by-NTHETA-by-NROT-by-NBETA complex*16 matrix
c
c  Our definition of complex spherical harmonics is
c
c  Ynm(theta,phi)= sqrt( 2n+1) sqrt((n-m)!/(n+m)!) 
c                  Pnm(cos theta) e^(im phi), 
c  Yn,-m(theta,phi) = sqrt( 2n+1) sqrt((n-m)!/(n+m)!) 
c                  Pnm(cos theta) e^(-im phi),   for m >= 0.
c       
c  Note that we do not include the Condon-Shortley phase (-1)^m, if m<0.
c
        complex *16 mpole(0:nterms,-nterms:nterms)

        real *8 beta(nbeta)
        real *8 rotmat(0:nterms,-nterms:nterms,-nterms:nterms,nbeta) 
        real *8 ctheta(ntheta),ynms(0:nterms,0:nterms,ntheta/2+1)
        complex *16 wsave(4*nphi+15)
c
        complex *16 grids(nphi,ntheta,nrot,nbeta)
c
        complex *16, allocatable :: marray(:,:) 
        complex *16, allocatable :: mptmp(:,:) 
c
        allocate( marray(0:nterms,-nterms:nterms) )
        allocate( mptmp(0:nterms,-nterms:nterms) )
c
        done=1
        pi=4*atan(done)
c
        do k=1,nbeta
        do i=1,nrot
        alpha=(i-1)*2*pi/dble(nrot)
        call mpolerotz_ifast_cmpl(nterms,mpole,alpha,mptmp)
        call rotviarecur3px_apply
     $     (beta(k),nterms,nterms,nterms,
     $     mptmp,nterms,marray,nterms,
     $     rotmat(0,-nterms,-nterms,k),nterms)

        call sphtrans_cmpl(
     $     nterms,marray(0,-nterms),
     $     nphi,ntheta,grids(1,1,i,k),
     $     ctheta,ynms,wsave)
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine rotgrid_fsr_real_init(nterms,nbeta,beta,rotmat)
        implicit real *8 (a-h,o-z)
c
c  Initialize rotgrid_fsr_real.
c
c  Fast, FFT-based algorithm for rotating real spherical harmonic grids
c  into pole locations (beta_j, alpha_k), 
c  where alpha_k = 2*pi * k/nrot, k=0..nrot-1, beta_j=1..nbeta
c
c  Input parameters:
c
c  nterms - the number of terms in spherical harmonics expansion
c  nbeta - the number of points in pole location meridian discretization
c  beta - angles for new pole locations beta_j, j=1..nbeta
c
c  Output parameters:
c
c  rotmat - The rotation operators for directions (beta_j,0). 
c
        real *8 beta(nbeta)
        real *8 rotmat(0:nterms,-nterms:nterms,-nterms:nterms,nbeta) 
        complex *16, allocatable :: rotmat_cmpl(:,:,:)
c
        allocate( rotmat_cmpl(0:nterms,-nterms:nterms,-nterms:nterms) )
c
        alpha=0
c
        do kk=1,nbeta
        call rotmat_projf90(beta(kk),alpha,nterms,rotmat_cmpl,nterms)
c
        do j=-nterms,nterms
        do i=-nterms,nterms
        do k=0,nterms
        rotmat(k,i,j,kk) = rotmat_cmpl(k,i,j)
        enddo
        enddo
        enddo

        enddo
c
        return
        end
c
c
c
c
c
        subroutine rotgrid_fsr_real(nterms,mpole,nphi,ntheta,
     $     nbeta,beta,nrot,grids,rotmat,ctheta,ynms,wsave)
        implicit real *8 (a-h,o-z)
c
c  Rotate the real spherical harmonic grids.
c
c  Fast, FFT-based algorithm for rotating real spherical harmonic grids
c  (dimensioned NPHI-by-NTHETA) into new pole locations (beta_j, alpha_k), 
c  where alpha_k = 2*pi * k/nrot, k=0..nrot-1,  beta_j, j=1..nbeta. 
c
c  GRIDS = rotgrid_fsr_real(NTERMS,MPOLE,NPHI,NTHETA,NBETA,BETA,NROT,...
c       ROTMAT,CTHETA,YNMS,WSAVE) rotates 
c  the real spherical harmonics expansion of degree NTERMS 
c  about the z-axis by degree ALPHA_K and about the y-axis by degree BETA.
c  into a collection of new pole locations (beta_j, alpha_k) in spherical 
c  coordinates (theta, phi), where alpha_k = 2*pi * k/nrot, k=0..nrot-1,
c  beta_j, j=1..nbeta.
c
c  The rotated poles form a uniformly spaced grid on lattitude \theta.
c
c  grids - function values on the rotated grids, 
c             NPHI-by-NTHETA-by-NROT-by-NBETA real*8 matrix
c
c      Input parameters:
c
c  nterms - the number of terms in spherical harmonics expansion
c  mpole - the coefficients of spherical harmonics expansion,
c                    complex*16 (0:nterms,0:nterms)
c  nphi - the number of points in latitude discretization (for spherical grid)
c  ntheta - the number of points in meridian discretization (for spherical grid)
c  ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
c  ynms - Legendre functions, real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
c  wsave - fftpack wsave array, complex *16 wsave(4*nphi+15)
c
c  nbeta - the number of points in pole location meridian discretization
c  beta - angles for new pole locations beta_j, j=1..nbeta
c  nrot - angles for new pole locations alpha_k = 2*pi * k/nrot, k=0..nrot-1.
c       
c  rotmat - The rotation operators for directions (beta_j,0). 
c 
c  [rotmat, ctheta, ynms, wsave] must be initialized via 
c              a preceding call to rotgrid_fsr_real_init
c
c      Output parameters:
c
c  grids - function values on the rotated grids, 
c             NPHI-by-NTHETA-by-NROT-by-NBETA real*8 matrix
c
c  Our definition of complex spherical harmonics is
c
c  Ynm(theta,phi)= sqrt( 2n+1) sqrt((n-m)!/(n+m)!) 
c                  Pnm(cos theta) e^(im phi), 
c  Yn,-m(theta,phi) = sqrt( 2n+1) sqrt((n-m)!/(n+m)!) 
c                  Pnm(cos theta) e^(-im phi),   for m >= 0.
c       
c  Note that we do not include the Condon-Shortley phase (-1)^m, if m<0.
c
        complex *16 mpole(0:nterms,0:nterms)

        real *8 beta(nbeta)
        real *8 rotmat(0:nterms,-nterms:nterms,-nterms:nterms,nbeta) 
        real *8 ctheta(ntheta),ynms(0:nterms,0:nterms,ntheta/2+1)
        complex *16 wsave(4*nphi+15)
c
        real *8 grids(nphi,ntheta,nrot,nbeta)
c
        complex *16, allocatable :: marray(:,:,:) 
c
        allocate( marray(0:nterms,0:nterms,nrot) )
c
        do k=1,nbeta
        call rot1lat_wfft_real
     $     (beta(k),nrot,nterms,
     $     nterms,nterms,mpole,nterms,marray,
     $     nterms,rotmat(0,-nterms,-nterms,k),nterms)
        do i=1,nrot
        call sphtrans_real(
     $     nterms,marray(0,0,i),
     $     nphi,ntheta,grids(1,1,i,k),
     $     ctheta,ynms,wsave)
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine rotgrid_dsr_real(nterms,mpole,nphi,ntheta,
     $     nbeta,beta,nrot,grids,rotmat,ctheta,ynms,wsave)
        implicit real *8 (a-h,o-z)
c
c  Rotate the real spherical harmonic grids.
c
c  Fast, direct algorithm for rotating real spherical harmonic grids
c  (dimensioned NPHI-by-NTHETA) into new pole locations (beta_j, alpha_k), 
c  where alpha_k = 2*pi * k/nrot, k=0..nrot-1,  beta_j, j=1..nbeta. 
c
c  GRIDS = rotgrid_fsr_real(NTERMS,MPOLE,NPHI,NTHETA,NBETA,BETA,NROT,...
c       ROTMAT,CTHETA,YNMS,WSAVE) rotates 
c  the real spherical harmonics expansion of degree NTERMS 
c  about the z-axis by degree ALPHA_K and about the y-axis by degree BETA.
c  into a collection of new pole locations (beta_j, alpha_k) in spherical 
c  coordinates (theta, phi), where alpha_k = 2*pi * k/nrot, k=0..nrot-1,
c  beta_j, j=1..nbeta.
c
c  The rotated poles form a uniformly spaced grid on lattitude \theta.
c
c  grids - function values on the rotated grids, 
c             NPHI-by-NTHETA-by-NROT-by-NBETA real*8 matrix
c
c      Input parameters:
c
c  nterms - the number of terms in spherical harmonics expansion
c  mpole - the coefficients of spherical harmonics expansion,
c                    complex*16 (0:nterms,0:nterms)
c  nphi - the number of points in latitude discretization (for spherical grid)
c  ntheta - the number of points in meridian discretization (for spherical grid)
c  ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
c  ynms - Legendre functions, real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
c  wsave - fftpack wsave array, complex *16 wsave(4*nphi+15)
c
c  nbeta - the number of points in pole location meridian discretization
c  beta - angles for new pole locations beta_j, j=1..nbeta
c  nrot - angles for new pole locations alpha_k = 2*pi * k/nrot, k=0..nrot-1.
c       
c  rotmat - The rotation operators for directions (beta_j,0). 
c 
c  [rotmat, ctheta, ynms, wsave] must be initialized via 
c              a preceding call to rotgrid_fsr_real_init
c
c      Output parameters:
c
c  grids - function values on the rotated grids, 
c             NPHI-by-NTHETA-by-NROT-by-NBETA real*8 matrix
c
c  Our definition of complex spherical harmonics is
c
c  Ynm(theta,phi)= sqrt( 2n+1) sqrt((n-m)!/(n+m)!) 
c                  Pnm(cos theta) e^(im phi), 
c  Yn,-m(theta,phi) = sqrt( 2n+1) sqrt((n-m)!/(n+m)!) 
c                  Pnm(cos theta) e^(-im phi),   for m >= 0.
c       
c  Note that we do not include the Condon-Shortley phase (-1)^m, if m<0.
c
        complex *16 mpole(0:nterms,0:nterms)

        real *8 beta(nbeta)
        real *8 rotmat(0:nterms,-nterms:nterms,-nterms:nterms,nbeta) 
        real *8 ctheta(ntheta),ynms(0:nterms,0:nterms,ntheta/2+1)
        complex *16 wsave(4*nphi+15)
c
        real *8 grids(nphi,ntheta,nrot,nbeta)
c
        complex *16, allocatable :: marray(:,:) 
        complex *16, allocatable :: mptmp(:,:) 
c
        allocate( marray(0:nterms,0:nterms) )
        allocate( mptmp(0:nterms,0:nterms) )
c
        done=1
        pi=4*atan(done)
c
        do k=1,nbeta
        do i=1,nrot

        alpha=(i-1)*2*pi/dble(nrot)
        call mpolerotz_ifast_real(nterms,mpole,alpha,mptmp)
        call rotviarecur3px_apply_real
     $     (beta(k),nterms,nterms,nterms,
     $     mptmp,nterms,marray,nterms,
     $     rotmat(0,-nterms,-nterms,k),nterms)

        call sphtrans_real(
     $     nterms,marray(0,0),
     $     nphi,ntheta,grids(1,1,i,k),
     $     ctheta,ynms,wsave)
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine rotgrid_vec(nphi,phi,ntheta,theta,fmodes,
     $     ngrid,nlat,xgrid,ygrid,zgrid,grids)
        implicit real *8 (a-h,o-z)
c
c       Fast rotation of the user-defined grid (xgrid,ygrid,zgrid).
c       This function contructs a set of grids, obtained by rotating 
c       the original grid by nphi uniformly spaced angles phi.
c       
c       Vectorized version, for multiple latitude grids.
c
c       The output is NGRID-by-NPHI-by-NLAT complex matrix, 
c       containing the rotated grids.
c       
c       Input parameters:
c       
c       nphi - the number of points in latitude discretization
c       phi - lattitude discretization angles
c       ntheta - the number of points in big circle discretization
c       theta - big circle discretization angles
c       fmodes - Fourier modes, NPHI-by-NTHETA complex matrix
c       ngrid - the number of points in the user-defined grid
c       nlat - the number latitudes/grids
c       xgrid - the x-coordinates of the grid, NGRID-BY-NLAT real matrix
c       ygrid - the y-coordinates of the grid, NGRID-BY-NLAT real matrix
c       zgrid - the z-coordinates of the grid, NGRID-BY-NLAT real matrix
c       
c       Output parameters:
c       
c       grids - function values at rotated grids, 
c                         NGRID-by-NPHI-by-NLAT complex matrix
c
c
        real *8 phi(nphi),theta(ntheta)
        complex *16 grids(ngrid,nphi,nlat)
        complex *16 fmodes(nphi,ntheta)
        real *8 xgrid(ngrid,nlat),ygrid(ngrid,nlat),zgrid(ngrid,nlat)

        do i=1,nlat
        call rotgrid(nphi,phi,ntheta,theta,fmodes,
     $     ngrid,xgrid(1,i),ygrid(1,i),zgrid(1,i),grids(1,1,i))
        enddo

        return
        end
c
c
c
c
c
        subroutine rotgrid_real_vec(nphi,phi,ntheta,theta,fmodes,
     $     ngrid,nlat,xgrid,ygrid,zgrid,grids)
        implicit real *8 (a-h,o-z)
c
c       Fast rotation of the user-defined grid (xgrid,ygrid,zgrid).
c       This function contructs a set of grids, obtained by rotating 
c       the original grid by nphi uniformly spaced angles phi.
c       
c       Vectorized version, for multiple latitude grids.
c
c       The output is NGRID-by-NPHI-by-NLAT real matrix, 
c       containing the rotated grids.
c       
c       Input parameters:
c       
c       nphi - the number of points in latitude discretization
c       phi - lattitude discretization angles
c       ntheta - the number of points in big circle discretization
c       theta - big circle discretization angles
c       fmodes - Fourier modes, NPHI-by-NTHETA complex matrix
c       ngrid - the number of points in the user-defined grid
c       nlat - the number latitudes/grids
c       xgrid - the x-coordinates of the grid, NGRID-BY-NLAT real matrix
c       ygrid - the y-coordinates of the grid, NGRID-BY-NLAT real matrix
c       zgrid - the z-coordinates of the grid, NGRID-BY-NLAT real matrix
c       
c       Output parameters:
c       
c       grids - function values at rotated grids, 
c                         NGRID-by-NPHI-by-NLAT real matrix
c
c
        real *8 phi(nphi),theta(ntheta)
        real *8 grids(ngrid,nphi,nlat)
        complex *16 fmodes(nphi,ntheta)
        real *8 xgrid(ngrid,nlat),ygrid(ngrid,nlat),zgrid(ngrid,nlat)

ccc        t1=second()
        do i=1,nlat
        call rotgrid_real(nphi,phi,ntheta,theta,fmodes,
     $     ngrid,xgrid(1,i),ygrid(1,i),zgrid(1,i),grids(1,1,i))
        enddo
ccc        t2=second()
c        call prini(6,13)
c        call prin2('in rotgrid_real_opt_vec, time=*', t2-t1, 1)
c        write(*,*) 'in rotgrid_real_opt_vec, time=*', t2-t1
        return
        end
c
c
c
c
c
        subroutine rotgrid_opt_vec(nphi,phi,ntheta,theta,fmodes,
     $     ngrid,nlat,xgrid,ygrid,zgrid,grids)
        implicit real *8 (a-h,o-z)
c
c       Fast rotation of the user-defined grid (xgrid,ygrid,zgrid).
c       This function contructs a set of grids, obtained by rotating 
c       the original grid by nphi uniformly spaced angles phi.
c       
c       Vectorized version, for multiple latitude grids.
c
c       The output is NGRID-by-NPHI-by-NLAT complex matrix, 
c       containing the rotated grids.
c       
c  Optimized for symmetric spherical grid, nphi is even, ngrid = nphi*ntheta
c  Note, this routine does not work for arbitrary grids, use rotgrid.
c
c       Input parameters:
c       
c       nphi - the number of points in latitude discretization
c       phi - lattitude discretization angles
c       ntheta - the number of points in big circle discretization
c       theta - big circle discretization angles
c       fmodes - Fourier modes, NPHI-by-NTHETA complex matrix
c       ngrid - the number of points in the user-defined grid
c       xgrid - the x-coordinates of the grid, NGRID-BY-NLAT real matrix
c       ygrid - the y-coordinates of the grid, NGRID-BY-NLAT real matrix
c       zgrid - the z-coordinates of the grid, NGRID-BY-NLAT real matrix
c       
c       Output parameters:
c       
c       grids - function values at rotated grids, 
c               NGRID-by-NPHI-by-NLAT complex matrix
c
c
        real *8 phi(nphi),theta(ntheta)
        complex *16 grids(ngrid,nphi,nlat)
        complex *16 fmodes(nphi,ntheta)
        real *8 xgrid(ngrid,nlat),ygrid(ngrid,nlat),zgrid(ngrid,nlat)

        do i=1,nlat
        call rotgrid_opt(nphi,phi,ntheta,theta,fmodes,
     $     ngrid,xgrid(1,i),ygrid(1,i),zgrid(1,i),grids(1,1,i))
        enddo

        return
        end
c
c
c
c
c
        subroutine rotgrid_real_opt_vec(nphi,phi,ntheta,theta,fmodes,
     $     ngrid,nlat,xgrid,ygrid,zgrid,grids)
        implicit real *8 (a-h,o-z)
c
c       Fast rotation of the user-defined grid (xgrid,ygrid,zgrid).
c       This function contructs a set of grids, obtained by rotating 
c       the original grid by nphi uniformly spaced angles phi.
c       
c       Vectorized version, for multiple latitude grids.
c
c       The output is NGRID-by-NPHI-by-NLAT real matrix, 
c       containing the rotated grids.
c       
c  Optimized for symmetric spherical grid, nphi is even, ngrid = nphi*ntheta
c  Note, this routine does not work for arbitrary grids, use rotgrid_real.
c
c       Input parameters:
c       
c       nphi - the number of points in latitude discretization
c       phi - lattitude discretization angles
c       ntheta - the number of points in big circle discretization
c       theta - big circle discretization angles
c       fmodes - Fourier modes, NPHI-by-NTHETA complex matrix
c       ngrid - the number of points in the user-defined grid
c       xgrid - the x-coordinates of the grid, NGRID-BY-NLAT real matrix
c       ygrid - the y-coordinates of the grid, NGRID-BY-NLAT real matrix
c       zgrid - the z-coordinates of the grid, NGRID-BY-NLAT real matrix
c       
c       Output parameters:
c       
c       grids - function values at rotated grids, 
c               NGRID-by-NPHI-by-NLAT real matrix
c
c
        real *8 phi(nphi),theta(ntheta)
        real *8 grids(ngrid,nphi,nlat)
        complex *16 fmodes(nphi,ntheta)
        real *8 xgrid(ngrid,nlat),ygrid(ngrid,nlat),zgrid(ngrid,nlat)

ccc        t1=second()
        do i=1,nlat
        call rotgrid_real_opt(nphi,phi,ntheta,theta,fmodes,
     $     ngrid,xgrid(1,i),ygrid(1,i),zgrid(1,i),grids(1,1,i))
        enddo
ccc        t2=second()
c        call prini(6,13)
c        call prin2('in rotgrid_real_opt_vec, time=*', t2-t1, 1)
c        write(*,*) 'in rotgrid_real_opt_vec, time=*', t2-t1
        return
        end
c
c
c
c
c
        subroutine rotnlat_fsr_cmpl_init(nterms,nbeta,beta,rotmat)
        implicit real *8 (a-h,o-z)
        real *8 beta(nbeta)
        real *8 rotmat(0:nterms,-nterms:nterms,-nterms:nterms,nbeta) 
        complex *16, allocatable :: rotmat_cmpl(:,:,:)

        allocate( rotmat_cmpl(0:nterms,-nterms:nterms,-nterms:nterms) )

        alpha=0
c
        do kk=1,nbeta
        call rotmat_projf90(beta(kk),alpha,nterms,rotmat_cmpl,nterms)
c
        do j=-nterms,nterms
        do i=-nterms,nterms
        do k=0,nterms
        rotmat(k,i,j,kk) = rotmat_cmpl(k,i,j)
        enddo
        enddo
        enddo

        enddo
c
        return
        end
c
c
c
c
c
        subroutine rotnlat_fsr_cmpl(nterms,mpole,
     $     nbeta,beta,nrot,marray,rotmat)
        implicit real *8 (a-h,o-z)
        complex *16 mpole(0:nterms,-nterms:nterms)

        real *8 beta(nbeta)
        real *8 rotmat(0:nterms,-nterms:nterms,-nterms:nterms,nbeta) 
c
        complex *16 marray(0:nterms,-nterms:nterms,nrot,nbeta)
c
        do k=1,nbeta
        call rot1lat_wfft
     $     (beta(k),nrot,nterms,
     $     nterms,nterms,mpole,nterms,marray(0,-nterms,1,k),
     $     nterms,rotmat(0,-nterms,-nterms,k),nterms)
        enddo
c
        return
        end
c
c
c
c
c
        subroutine mpolerotz_ifast_cmpl(nterms,mpole,phi,mpout)
        implicit real *8 (a-h,o-z)
c
c       This subroutine performs the rotation of the multipole expansion
c       about the z-axis
c
c     Input parameters:
c
c     nterms  (integer *4)  : number of terms in the multipole expansion
c     mpole  (complex *16)  : the multipole expansion
c     phi        (real *8)  : the angle of rotation
c
c     Output parameters:
c 
c     mpout  (complex *16)  : the rotated multipole expansion
c
c
        complex *16 mpole(0:nterms,-nterms:nterms)
        complex *16 mpout(0:nterms,-nterms:nterms)
        complex *16 ephi,ima,cd
        data ima/(0.0d0,1.0d0)/
c
c
c----- a rotation of PHI radians about the Z-axis 
c
c
        cd=exp(+ima*phi)
        ephi=1

        do m=0,nterms
        do l=m,nterms
        mpout(l,+m)=mpole(l,+m)*ephi
        mpout(l,-m)=mpole(l,-m)*dconjg(ephi)
        enddo
        ephi=ephi*cd
        enddo
c
        return
        end
c
c
c
c
c
        subroutine mpolerotz_ifast_real(nterms,mpole,phi,mpout)
        implicit real *8 (a-h,o-z)
c
c       This subroutine performs the rotation of the multipole expansion
c       about the z-axis
c
c     Input parameters:
c
c     nterms  (integer *4)  : number of terms in the multipole expansion
c     mpole  (complex *16)  : the multipole expansion
c     phi        (real *8)  : the angle of rotation
c
c     Output parameters:
c 
c     mpout  (complex *16)  : the rotated multipole expansion
c
c
        complex *16 mpole(0:nterms,0:nterms)
        complex *16 mpout(0:nterms,0:nterms)
        complex *16 ephi,ima,cd
        data ima/(0.0d0,1.0d0)/
c
c
c----- a rotation of PHI radians about the Z-axis 
c
c
        cd=exp(+ima*phi)
        ephi=1

        do m=0,nterms
        do l=m,nterms
        mpout(l,+m)=mpole(l,+m)*ephi
        enddo
        ephi=ephi*cd
        enddo
c
        return
        end
