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
c     ROTATION OF SPHERICAL GRIDS  (FORTRAN 77 AND 90 VERSIONS).
c


        subroutine init_grid_arb_single
     $     (nphi,ntheta, phi,theta,xs,ys,zs,ws)
        implicit real *8 (a-h,o-z)
c
c       Precompute parameters for constructing a spherical grid
c       NOTE: we are using meridians theta = [0..pi] in this function
c
c       Arbitrary grid
c
        real *8 phi(nphi),theta(ntheta)
        real *8 xs(nphi),ys(nphi),zs(ntheta),ws(ntheta)

        done=1
        pi=4*atan(done)

        do i=1,nphi
        phi(i)=(i-1)*2*pi/dble(nphi)
        xs(i)=cos(phi(i))
        ys(i)=sin(phi(i))
        enddo

        do i=1,ntheta
        theta(i)=(i-0.5d0)*1*pi/dble(ntheta)
        zs(i)=-cos(theta(i))
        ws(i)=2/dble(ntheta)
        if( theta(i) .gt. pi ) theta(i)=theta(i)-2*pi
        enddo

        return
        end
c
c
c
c
c
        subroutine init_grid_cheb_single
     $     (nphi,ntheta, phi,theta,xs,ys,zs,ws)
        implicit real *8 (a-h,o-z)
c
c       Precompute parameters for constructing a spherical grid
c       NOTE: we are using meridians theta = [0..pi] in this function
c
c       Chebychev grid
c
        real *8 phi(nphi),theta(ntheta)
        real *8 xs(nphi),ys(nphi),zs(ntheta),ws(ntheta)

        done=1
        pi=4*atan(done)

        do i=1,nphi
        phi(i)=(i-1)*2*pi/dble(nphi)
        xs(i)=cos(phi(i))
        ys(i)=sin(phi(i))
        enddo

        do i=1,ntheta
        theta(i)=(i-0.5d0)*1*pi/dble(ntheta)
        zs(i)=-cos(theta(i))
        ws(i)=2/dble(ntheta)
        if( theta(i) .gt. pi ) theta(i)=theta(i)-2*pi
        enddo

        return
        end
c
c
c
c
c
        subroutine init_grid_cheb_double
     $     (nphi,ntheta, phi,theta,xs,ys,zs,ws)
        implicit real *8 (a-h,o-z)
c
c       Precompute parameters for constructing a spherical grid
c       NOTE: we are using big circles theta = [-pi..pi] in this function
c
c       Chebychev grid
c
        real *8 phi(nphi),theta(ntheta)
        real *8 xs(nphi),ys(nphi),zs(ntheta),ws(ntheta)

        done=1
        pi=4*atan(done)

        do i=1,nphi
        phi(i)=(i-1)*2*pi/dble(nphi)
        xs(i)=cos(phi(i))
        ys(i)=sin(phi(i))
        enddo

        do i=1,ntheta
        theta(i)=(i-0.5d0)*2*pi/dble(ntheta)
        zs(i)=-cos(theta(i))
        ws(i)=4/dble(ntheta)
        if( theta(i) .gt. pi ) theta(i)=theta(i)-2*pi
        enddo

        return
        end
c
c
c
c
c
        subroutine init_grid_lege_single
     $     (nphi,ntheta, phi,theta,xs,ys,zs,ws)
        implicit real *8 (a-h,o-z)
c
c       Precompute parameters for constructing a spherical grid
c       NOTE: we are using meridians theta = [0..pi] in this function
c
c       Legendre grid
c
        real *8 phi(nphi),theta(ntheta)
        real *8 xs(nphi),ys(nphi),zs(ntheta),ws(ntheta)

        done=1
        pi=4*atan(done)

        do i=1,nphi
        phi(i)=(i-1)*2*pi/dble(nphi)
        xs(i)=cos(phi(i))
        ys(i)=sin(phi(i))
        enddo

        ifwhts=1
        call legewhts(ntheta,zs,ws,ifwhts)

        do i=1,ntheta
        zs(i)=-zs(i)
        enddo

        do i=1,ntheta
        theta(i)=acos(zs(i))
        if( theta(i) .gt. pi ) theta(i)=theta(i)-2*pi
        enddo

        return
        end
c
c
c
c
c
        subroutine init_grid_lege_double
     $     (nphi,ntheta, phi,theta,xs,ys,zs,ws)
        implicit real *8 (a-h,o-z)
c
c       Precompute parameters for constructing a spherical grid
c       NOTE: we are using big circles theta = [-pi..pi] in this function
c
c       Legendre grid
c
        real *8 phi(nphi),theta(ntheta)
        real *8 xs(nphi),ys(nphi),zs(ntheta),ws(ntheta)

        done=1
        pi=4*atan(done)

        do i=1,nphi
        phi(i)=(i-1)*2*pi/dble(nphi)
        xs(i)=cos(phi(i))
        ys(i)=sin(phi(i))
        enddo

        ifwhts=1
        call legewhts(ntheta/2,zs,ws,ifwhts)

        do i=1,ntheta/2
        zs(i)=-zs(i)
        enddo

        do i=1,ntheta/2
        zs(ntheta-i+1)=zs(i)
        ws(ntheta-i+1)=ws(i)
        enddo

        do i=1,ntheta/2
        theta(i)=acos(zs(i))
        theta(ntheta-i+1)=2*pi-theta(i)
        enddo

        do i=1,ntheta
        if( theta(i) .gt. pi ) theta(i)=theta(i)-2*pi
        enddo

        return
        end
c
c
c
c
c
        subroutine xyz_grid(beta,nphi,xs,ys,ntheta,zs,theta, 
     $     xgrid,ygrid,zgrid)
        implicit real *8 (a-h,o-z)
c
c       Construct a spherical grid and rotate it by angle beta around y-axis
c
c       Input parameters:
c
c       beta - the angle of rotation around y-axis
c       nphi, ntheta - the number of discretization points in phi and theta
c       phi,theta,xs,ys,zs - must be constructed via a preceding call to 
c                 either init_grid_single or init_grid_double
c
c       Output parameters:
c       
c       xgrid - the x-coordinates of the constructed grid
c       ygrid - the y-coordinates of the constructed grid
c       zgrid - the z-coordinates of the constructed grid
c       
c
        real *8 phi(nphi),theta(ntheta)
        real *8 xs(nphi),ys(nphi),zs(ntheta)
        real *8 xgrid(nphi,ntheta)
        real *8 ygrid(nphi,ntheta)
        real *8 zgrid(nphi,ntheta)

        cosbeta=cos(beta)
        sinbeta=sin(beta)

        do j=1,ntheta
        costheta=cos(theta(j))
        sintheta=sin(theta(j))
        do i=1,nphi
        z=costheta
        x=xs(i)*sintheta
        y=ys(i)*sintheta
c  counter-clockwise rotation: positive angle
c        xgrid(i,j) = x*cosbeta - z*sinbeta
c        ygrid(i,j) = y
c        zgrid(i,j) = x*sinbeta + z*cosbeta
c  clockwise rotation: negative angle
        xgrid(i,j) = x*cosbeta + z*sinbeta
        ygrid(i,j) = y
        zgrid(i,j) =-x*sinbeta + z*cosbeta
        enddo
        enddo

        return
        end
c
c
c
c
c
        subroutine xyz_grid_a(beta,nphi,xs,ys,ntheta,zs,theta,alpha,
     $     xgrid,ygrid,zgrid)
        implicit real *8 (a-h,o-z)
c
c       Construct a spherical grid and rotate it by angle beta around y-axis
c       then, rotate it by angle alpha around z-axis
c
c       Input parameters:
c
c       beta - the angle of rotation around y-axis
c       alpha - the angle of rotation around z-axis
c       nphi, ntheta - the number of discretization points in phi and theta
c       phi,theta,xs,ys,zs - must be constructed via a preceding call to 
c                 either init_grid_single or init_grid_double
c
c       Output parameters:
c       
c       xgrid - the x-coordinates of the constructed grid
c       ygrid - the y-coordinates of the constructed grid
c       zgrid - the z-coordinates of the constructed grid
c       
c
        real *8 phi(nphi),theta(ntheta)
        real *8 xs(nphi),ys(nphi),zs(ntheta)
        real *8 xgrid(nphi,ntheta)
        real *8 ygrid(nphi,ntheta)
        real *8 zgrid(nphi,ntheta)

        cosbeta=cos(beta)
        sinbeta=sin(beta)

        cosalpha=cos(alpha)
        sinalpha=sin(alpha)

        do j=1,ntheta
        costheta=cos(theta(j))
        sintheta=sin(theta(j))
        do i=1,nphi
        z=costheta
        x=xs(i)*sintheta
        y=ys(i)*sintheta
c  counter-clockwise rotation: positive angle
c        x1 = x*cosbeta - z*sinbeta
c        y1 = y
c        z1 = x*sinbeta + z*cosbeta
c  clockwise rotation: negative angle
        x1 = x*cosbeta + z*sinbeta
        y1 = y
        z1 =-x*sinbeta + z*cosbeta
        xgrid(i,j) = x1*cosalpha - y1*sinalpha
        ygrid(i,j) = x1*sinalpha + y1*cosalpha
        zgrid(i,j) = z1
        enddo
        enddo

        return
        end
c
c
c
c
c
        subroutine rotgrid3phi(x,y,phi)
        implicit real *8 (a-h,o-z)
c
        if( abs(x) .eq. 0 .and. abs(y) .eq. 0 ) then
        phi = 0
        else
        phi = atan2(y,x)
        endif
c
        return
        end
c
c
c
c
c
        subroutine rotgrid3phitheta(x,y,z,phi,theta)
        implicit real *8 (a-h,o-z)
c
        if( abs(x) .eq. 0 .and. abs(y) .eq. 0 ) then
        phi = 0
        else
        phi = atan2(y,x)
        endif
c
        r=sqrt(x**2+y**2)
        theta = atan2(r,z)
c
        return
        end
c
c
c
c
c
        subroutine rotgrid3theta(z,theta)
        implicit real *8 (a-h,o-z)
c
        if( z .ge. +1 ) then
        theta = 0
        else if( z .le. -1 ) then
        done = 1
        pi = 4*atan(done)
        theta = pi
        else
        theta = acos(z)
        endif
c
        return
        end
c
c
c
c
c
        subroutine rotgridi(nphi,phi,ntheta,theta,fgrid,fmodes)
        implicit real *8 (a-h,o-z)
c       
c       Get Fourier modes for a user-defined function on a spherical grid.
c       
c       The output is NPHI-by-NTHETA complex matrix, which contains
c       the values of the Fourier coefficients. Note, that the user data
c       must be specified as NPHI-by-NTHETA complex matrix.
c
c
c       Input parameters:
c       
c       nphi - the number of points in lattitude discretization
c       phi - lattitude discretization angles
c       ntheta - the number of points in big circle discretization
c       theta - big circle discretization angles
c       fgrid - user data, NPHI-by-NTHETA complex matrix
c
c       Output parameters:
c
c       fmodes - Fourier modes, NPHI-by-NTHETA complex matrix
c
c
        real *8 phi(nphi),theta(ntheta)
        complex *16 fgrid(nphi,ntheta),fmodes(nphi,ntheta)
        real *8, allocatable :: phi_grid(:,:)
        real *8, allocatable :: theta_grid(:,:)
        complex *16, allocatable :: wsave(:)
        real *8, allocatable :: phase(:)

        eps=1d-12 * 1e0

c
c Optionally, use nuFFT in R^2 to compare timings
c
        if_use_nufft2d = 1
        if( if_use_nufft2d .eq. 1 ) then
c
c Get Fourier coefficients for lattitudes
c Get Fourier coefficients for big circles
c
        allocate( phi_grid(nphi,ntheta) )
        allocate( theta_grid(nphi,ntheta) )

        do j=1,ntheta
        do i=1,nphi
        phi_grid(i,j)=phi(i)
        theta_grid(i,j)=theta(j)
        enddo
        enddo

        iflag=-1
        call nufft2d1f90(nphi*ntheta,phi_grid,theta_grid,fgrid,
     $     iflag,eps,nphi,ntheta,fmodes,ier)

ccc        call prin2('fmodes=*',fmodes,2*nphi*ntheta)
c        do i=1,ntheta
c        call prin2('fmodes=*',fmodes(1,i),2*nphi)
c        enddo

ccc     test symmetries for real valued data
c        write(*,*) nphi, 5, nphi-3, nphi/2
c        write(*,*) (fmodes(5,k), k=1,ntheta)
c        write(*,*)
c        write(*,*) (fmodes(nphi-3,k), k=1,ntheta)

        endif

        return
        end
c
c
c
c
c
        subroutine rotgrid(nphi,phi,ntheta,theta,fmodes,
     $     ngrid,xgrid,ygrid,zgrid,grids)
        implicit real *8 (a-h,o-z)
c
c       Fast rotation of the user-defined grid (xgrid,ygrid,zgrid).
c       This function contructs a set of grids, obtained by rotating 
c       the original grid by nphi uniformly spaced angles phi.
c       
c       The output is NGRID-by-NPHI complex matrix, 
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
c       xgrid - the x-coordinates of the grid
c       ygrid - the y-coordinates of the grid
c       zgrid - the z-coordinates of the grid
c       
c       Output parameters:
c       
c       grids - function values at rotated grids, NGRID-by-NPHI complex matrix
c
        real *8 phi(nphi),theta(ntheta)
        complex *16 grids(ngrid,nphi),fmodes(nphi,ntheta)
        real *8 xgrid(ngrid),ygrid(ngrid),zgrid(ngrid)
        real *8, allocatable :: zang(:)
        real *8, allocatable :: phase(:)
        real *8, allocatable :: phi_grid(:,:)
        real *8, allocatable :: theta_grid(:,:)

        complex *16, allocatable :: fmodes_tmp(:,:)
        complex *16, allocatable :: grid_tmp(:)
        complex *16, allocatable :: wsave(:)

        eps=1d-12 * 1e0

c
c ...and the phases and angles to be interpolated are
c
        allocate(zang(ngrid))
        allocate(phase(ngrid))
c
        do i=1,ngrid        
        call rotgrid3phitheta(xgrid(i),ygrid(i),zgrid(i),
     $     phase(i),zang(i))
        enddo
c
ccc        call prin2('zang=*',zang,ngrid)
ccc        call prin2('phase=*',phase,ngrid)


c
c Optionally, use nuFFT in R^2 to compare timings
c
        if_use_nufft2d = 0
        if( if_use_nufft2d .eq. 1 ) then
c
c Call nonuniform FFT in R^2 to get all grids simultaneously
c

        allocate( phi_grid(ngrid,nphi) )
        allocate( theta_grid(ngrid,nphi) )

        do i=1,ngrid
        do j=1,nphi
        phi_grid(i,j)=phi(j)+phase(i)
        theta_grid(i,j)=zang(i)
        enddo
        enddo

        iflag=+1
        call nufft2d2f90(ngrid*nphi,phi_grid,theta_grid,grids,
     $     iflag,eps,nphi,ntheta,fmodes,ier)

        return
        endif


        if_use_nufft1d = 1
        if( if_use_nufft1d .eq. 1 ) then
c
c Hybrid method
c Use nuFFT 1d in theta, and nuFFT/FFT 1d in phi
c

c
c Get interpolated values in theta
        iflag=+1
        allocate( fmodes_tmp(ntheta,nphi) )

        if( 1 .eq. 2 ) then
c
c Reference nuFFT in R^1
c        do k=1,nphi
c        do i=1,ntheta
c        fmodes_tmp(i,k)=fmodes(k,i)
c        enddo
c        call nufft1d2f90(ngrid,zang,grids(1,k),
c     $     iflag,eps,ntheta,fmodes_tmp(1,k),ier)
c        enddo

        endif

        if( 2 .eq. 2 ) then
c
c Vectorized nuFFT in R^1
        do k=1,nphi
        do i=1,ntheta
        fmodes_tmp(i,k)=fmodes(k,i)
        enddo
        enddo

c        do i=1,nphi
c        call prinf('i=*',i,1)
c        call prin2('fmodes_tmp=*',fmodes_tmp(1,i),2*ntheta)
c        enddo
c        stop

        call nufft1d2vf90(nphi,ngrid,zang,grids,
     $     iflag,eps,ntheta,fmodes_tmp,ier)

c        call prin2('zang=*',zang,ngrid)
c        do i=1,nphi
c        call prin2('grids=*',grids(1,i),2*ngrid)
c        enddo
c        stop

        endif

c
c Shift the grids on latitudes (apply phase in Fourier domain)
        call shiftphase(grids,ngrid,nphi,phase)

c
c FFT data on latitudes, apply inverse fftshift, even nphi only
c octave:9> ifftshift( 1:5 )
c ans =
c
c    3   4   5   1   2
c
c octave:10> ifftshift( 1:6 )
c ans =
c
c    4   5   6   1   2   3

        allocate( grid_tmp(nphi) )
        allocate( wsave(4*nphi+15) )
        call zffti(nphi,wsave)

        do i=1,ngrid
        do j=1,nphi
        j0=j+nphi/2
        if( j0 .gt. nphi ) j0=j0-nphi
        grid_tmp(j0)=grids(i,j)
        enddo
        call zfftb(nphi,grid_tmp,wsave)
        do j=1,nphi
        grids(i,j)=grid_tmp(j)
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
        subroutine rotgrid_opt(nphi,phi,ntheta,theta,fmodes,
     $     ngrid,xgrid,ygrid,zgrid,grids)
        implicit real *8 (a-h,o-z)
c
c       Fast rotation of the user-defined grid (xgrid,ygrid,zgrid).
c       This function contructs a set of grids, obtained by rotating 
c       the original grid by nphi uniformly spaced angles phi.
c       
c       Assume symmetric spherical grid, nphi is even, ngrid = nphi*ntheta
c
c       The output is NGRID-by-NPHI complex matrix, 
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
c       xgrid - the x-coordinates of the grid
c       ygrid - the y-coordinates of the grid
c       zgrid - the z-coordinates of the grid
c       
c       Output parameters:
c       
c       grids - function values at rotated grids, NGRID-by-NPHI complex matrix
c
        real *8 phi(nphi),theta(ntheta)
        complex *16 grids(ngrid,nphi),fmodes(nphi,ntheta)
        real *8 xgrid(ngrid),ygrid(ngrid),zgrid(ngrid)
        real *8, allocatable :: zang(:)
        real *8, allocatable :: phase(:)
        real *8, allocatable :: phi_grid(:,:)
        real *8, allocatable :: theta_grid(:,:)

        complex *16, allocatable :: fmodes_tmp(:,:)
        complex *16, allocatable :: grid_tmp(:)
        complex *16, allocatable :: wsave(:)

        real *8, allocatable :: zang_tmp(:,:)
        complex *16, allocatable :: grids_tmp(:,:,:)

        eps=1d-12 * 1e0

c
c ...and the phases and angles to be interpolated are
c
        allocate(zang(ngrid))
        allocate(phase(ngrid))
c
        do i=1,ngrid        
        call rotgrid3phitheta(xgrid(i),ygrid(i),zgrid(i),
     $     phase(i),zang(i))
        enddo
c
ccc        call prin2('zang=*',zang,ngrid)
ccc        call prin2('phase=*',phase,ngrid)


c
c Optionally, use nuFFT in R^2 to compare timings
c
        if_use_nufft2d = 0
        if( if_use_nufft2d .eq. 1 ) then
c
c Call nonuniform FFT in R^2 to get all grids simultaneously
c

        allocate( phi_grid(ngrid,nphi) )
        allocate( theta_grid(ngrid,nphi) )

        do i=1,ngrid
        do j=1,nphi
        phi_grid(i,j)=phi(j)+phase(i)
        theta_grid(i,j)=zang(i)
        enddo
        enddo

        iflag=+1
        call nufft2d2f90(ngrid*nphi,phi_grid,theta_grid,grids,
     $     iflag,eps,nphi,ntheta,fmodes,ier)

        return
        endif


        if_use_nufft1d = 1
        if( if_use_nufft1d .eq. 1 ) then
c
c Hybrid method
c Use nuFFT 1d in theta, and nuFFT/FFT 1d in phi
c

c
c Get interpolated values in theta
        iflag=+1
        allocate( fmodes_tmp(ntheta,nphi) )

        if( 1 .eq. 2 ) then
c
c Reference nuFFT in R^1
c        do k=1,nphi
c        do i=1,ntheta
c        fmodes_tmp(i,k)=fmodes(k,i)
c        enddo
c        call nufft1d2f90(ngrid,zang,grids(1,k),
c     $     iflag,eps,ntheta,fmodes_tmp(1,k),ier)
c        enddo

        endif

        if( 1 .eq. 2 ) then
c
c Vectorized nuFFT in R^1
        do k=1,nphi
        do i=1,ntheta
        fmodes_tmp(i,k)=fmodes(k,i)
        enddo
        enddo

c        do i=1,nphi
c        call prinf('i=*',i,1)
c        call prin2('fmodes_tmp=*',fmodes_tmp(1,i),2*ntheta)
c        enddo
c        stop

        call nufft1d2vf90(nphi,ngrid,zang,grids,
     $     iflag,eps,ntheta,fmodes_tmp,ier)

c        call prin2('zang=*',zang,ngrid)
c        do i=1,nphi
c        call prin2('grids=*',grids(1,i),2*ngrid)
c        enddo
c        stop

        endif

        if( 2 .eq. 2 ) then
c
c Assume symmetric grid, nphi is even        
c Vectorized nuFFT in R^1
        do k=1,nphi
        do i=1,ntheta
        fmodes_tmp(i,k)=fmodes(k,i)
        enddo
        enddo

        allocate(zang_tmp(nphi/2+1,ntheta/2))
        allocate(grids_tmp(nphi/2+1,ntheta/2,nphi))
        kk=1
        do i=1,ntheta/2
        do k=1,nphi/2+1
        zang_tmp(k,i)=zang(kk)
        kk=kk+1
        enddo
        kk=1+i*nphi
        enddo

c        do i=1,nphi/2
c        call prinf('i=*',i,1)
c        call prin2('fmodes_tmp=*',fmodes_tmp(1,i),2*ntheta)
c        enddo
c        stop

        ngrid_tmp = (nphi/2+1)*ntheta/2
        call nufft1d2vf90(nphi,ngrid_tmp,zang_tmp,grids_tmp,
     $     iflag,eps,ntheta,fmodes_tmp,ier)


c        write(*,*) ngrid_tmp, ntheta, nphi
c        do i=1,nphi
c        call prinf('i=*',i,1)
c        call prin2('grids_tmp=*',grids_tmp(1,1,i),2*ngrid_tmp)
c        enddo

        do j=1,nphi
        kk=1
        do i=1,ntheta/2
        do k=1,nphi/2+1
        grids(kk,j)=grids_tmp(k,i,j)
        kk=kk+1
        enddo
        do k=nphi/2+2,nphi
        grids(kk,j)=grids_tmp(nphi+2-k,i,j)
        kk=kk+1
        enddo
        enddo
        enddo

c        write(*,*) ngrid, ntheta, nphi
c        do i=1,nphi
c        call prinf('i=*',i,1)
c        call prin2('grids=*',grids(1,i),2*ngrid)
c        enddo
c        stop

        endif

c
c Shift the grids on latitudes (apply phase in Fourier domain)
        call shiftphase(grids,ngrid,nphi,phase)

c
c FFT data on latitudes, apply inverse fftshift, even nphi only
c octave:9> ifftshift( 1:5 )
c ans =
c
c    3   4   5   1   2
c
c octave:10> ifftshift( 1:6 )
c ans =
c
c    4   5   6   1   2   3

        allocate( grid_tmp(nphi) )
        allocate( wsave(4*nphi+15) )
        call zffti(nphi,wsave)

        do i=1,ngrid
        do j=1,nphi
        j0=j+nphi/2
        if( j0 .gt. nphi ) j0=j0-nphi
        grid_tmp(j0)=grids(i,j)
        enddo
        call zfftb(nphi,grid_tmp,wsave)
        do j=1,nphi
        grids(i,j)=grid_tmp(j)
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
        subroutine rotgridi_real(nphi,phi,ntheta,theta,fgrid,fmodes)
        implicit real *8 (a-h,o-z)
c       
c       Get Fourier modes for a user-defined function on a spherical grid.
c       
c       The output is NPHI-by-NTHETA complex matrix, which contains
c       the values of the Fourier coefficients. Note, that the user data
c       must be specified as NPHI-by-NTHETA complex matrix.
c
c
c       Input parameters:
c       
c       nphi - the number of points in lattitude discretization
c       phi - lattitude discretization angles
c       ntheta - the number of points in big circle discretization
c       theta - big circle discretization angles
c       fgrid - user data, NPHI-by-NTHETA complex matrix
c
c       Output parameters:
c
c       fmodes - Fourier modes, NPHI-by-NTHETA complex matrix
c
c
        real *8 phi(nphi),theta(ntheta)
        real *8 fgrid(nphi,ntheta)
        complex *16 fmodes(nphi,ntheta)
        real *8, allocatable :: phi_grid(:,:)
        real *8, allocatable :: theta_grid(:,:)
        complex *16, allocatable :: fgrid_cmpl(:,:)

        eps=1d-12 * 1e0

c
c Optionally, use nuFFT in R^2 to compare timings
c
        if_use_nufft2d = 1
        if( if_use_nufft2d .eq. 1 ) then
c
c Get Fourier coefficients for lattitudes
c Get Fourier coefficients for big circles
c
        allocate( phi_grid(nphi,ntheta) )
        allocate( theta_grid(nphi,ntheta) )

        allocate( fgrid_cmpl(nphi,ntheta) )

        do j=1,ntheta
        do i=1,nphi
        phi_grid(i,j)=phi(i)
        theta_grid(i,j)=theta(j)
        fgrid_cmpl(i,j)=fgrid(i,j)
        enddo
        enddo

        iflag=-1
        call nufft2d1f90(nphi*ntheta,phi_grid,theta_grid,fgrid_cmpl,
     $     iflag,eps,nphi,ntheta,fmodes,ier)

ccc     test symmetries for real valued data
c        write(*,*) nphi, 5, nphi-3, nphi/2
c        write(*,*) (fmodes(5,k), k=1,ntheta)
c        pause
c        write(*,*)
c        write(*,*) (fmodes(nphi-3,k), k=1,ntheta)
c        pause
        endif

        return
        end
c
c
c
c
c
        subroutine rotgrid_real(nphi,phi,ntheta,theta,fmodes,
     $     ngrid,xgrid,ygrid,zgrid,grids)
        implicit real *8 (a-h,o-z)
c
c       Fast rotation of the user-defined grid (xgrid,ygrid,zgrid).
c       This function contructs a set of grids, obtained by rotating 
c       the original grid by nphi uniformly spaced angles phi.
c       
c       The output is NGRID-by-NPHI real matrix, 
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
c       xgrid - the x-coordinates of the grid
c       ygrid - the y-coordinates of the grid
c       zgrid - the z-coordinates of the grid
c       
c       Output parameters:
c       
c       grids - function values at rotated grids, NGRID-by-NPHI complex matrix
c
        real *8 phi(nphi),theta(ntheta)
        real *8 grids(ngrid,nphi)
        complex *16 fmodes(nphi,ntheta)
        real *8 xgrid(ngrid),ygrid(ngrid),zgrid(ngrid)
        real *8, allocatable :: zang(:)
        real *8, allocatable :: phase(:)
        real *8, allocatable :: phi_grid(:,:)
        real *8, allocatable :: theta_grid(:,:)

        complex *16, allocatable :: fmodes_tmp(:,:)
        complex *16, allocatable :: grid_ctmp(:)
        real *8, allocatable :: grid_rtmp(:)
        complex *16, allocatable :: wsave(:)

        complex *16, allocatable :: grids_cmpl(:,:)

        eps=1d-12 * 1e0

c
c ...and the phases and angles to be interpolated are
c
        allocate(zang(ngrid))
        allocate(phase(ngrid))
c
        do i=1,ngrid        
        call rotgrid3phitheta(xgrid(i),ygrid(i),zgrid(i),
     $     phase(i),zang(i))
        enddo
c
ccc        call prin2('zang=*',zang,ngrid)
ccc        call prin2('phase=*',phase,ngrid)


c
c Optionally, use nuFFT in R^2 to compare timings
c
        if_use_nufft2d = 0
        if( if_use_nufft2d .eq. 1 ) then
c
c Call nonuniform FFT in R^2 to get all grids simultaneously
c

        allocate( phi_grid(ngrid,nphi) )
        allocate( theta_grid(ngrid,nphi) )

        allocate( grids_cmpl(ngrid,nphi) )

        do j=1,nphi
        do i=1,ngrid
        phi_grid(i,j)=phi(j)+phase(i)
        theta_grid(i,j)=zang(i)
        enddo
        enddo

        iflag=+1
        call nufft2d2f90(ngrid*nphi,phi_grid,theta_grid,grids_cmpl,
     $     iflag,eps,nphi,ntheta,fmodes,ier)

        do j=1,nphi
        do i=1,ngrid
        grids(i,j)=grids_cmpl(i,j)
        enddo
        enddo

        return
        endif


        if_use_nufft1d = 1
        if( if_use_nufft1d .eq. 1 ) then
c
c Hybrid method
c Use nuFFT 1d in theta, and nuFFT/FFT 1d in phi
c

c
c Get interpolated values in theta
        iflag=+1
        allocate( fmodes_tmp(ntheta,nphi) )

        allocate( grids_cmpl(ngrid,nphi) )

        if( 1 .eq. 2 ) then
c
c Reference nuFFT in R^1
c        do k=1,nphi
c        do i=1,ntheta
c        fmodes_tmp(i,k)=fmodes(k,i)
c        enddo
c        call nufft1d2f90(ngrid,zang,grids_cmpl(1,k),
c     $     iflag,eps,ntheta,fmodes_tmp(1,k),ier)
c        enddo
c        
        endif

        if( 1 .eq. 2 ) then
c
c Vectorized nuFFT in R^1
c        do k=1,nphi
c        do i=1,ntheta
c        fmodes_tmp(i,k)=fmodes(k,i)
c        enddo
c        enddo
c        call nufft1d2vf90(nphi,ngrid,zang,grids_cmpl,
c     $     iflag,eps,ntheta,fmodes_tmp,ier)

c        do i=1,nphi
c        call prinf('i=*',i,1)
c        call prin2('grids_cmpl=*',grids_cmpl(1,i),2*ngrid)
c        enddo
c        stop

        endif

        if( 2 .eq. 2 ) then
c
c Vectorized nuFFT in R^1, with real symmetries
        do k=1,nphi/2+1
        do i=1,ntheta
        fmodes_tmp(i,k)=fmodes(k,i)
        enddo
        enddo

c        do i=1,nphi/2
c        call prinf('i=*',i,1)
c        call prin2('fmodes_tmp=*',fmodes_tmp(1,i),2*ntheta)
c        enddo
c        stop

        call nufft1d2vf90(nphi/2+1,ngrid,zang,grids_cmpl,
     $     iflag,eps,ntheta,fmodes_tmp,ier)


c        do i=1,nphi/2+1
c        call prinf('i=*',i,1)
c        call prin2('grids_cmpl=*',grids_cmpl(1,i),2*ngrid)
c        enddo
c        stop

        endif

        if( 1 .eq. 2 ) then
c
c initialize the fourier grid for shiftphase and shiftphase_real routine
c
        do i=nphi/2+2,nphi
cc        write(*,*) i,nphi-i+2
        i0=nphi-i+2
        do j=1,ngrid
        grids_cmpl(j,i)=dconjg(grids_cmpl(j,i0))
        enddo
        enddo

        endif

c
c Shift the grids on latitudes (apply phase in Fourier domain)
ccc        call shiftphase(grids_cmpl,ngrid,nphi,phase)
ccc        call shiftphase_real(grids_cmpl,ngrid,nphi,phase)
        call shiftphase_real_alt(grids_cmpl,ngrid,nphi,phase)


c
c FFT data on latitudes, apply inverse fftshift, even nphi only
c octave:9> ifftshift( 1:5 )
c ans =
c
c    3   4   5   1   2
c
c octave:10> ifftshift( 1:6 )
c ans =
c
c    4   5   6   1   2   3

        allocate( grid_ctmp(nphi) )
        allocate( grid_rtmp(nphi) )
        allocate( wsave(4*nphi+15) )


        if ( 1 .eq. 2 ) then
c
c       ... complex version, needs shiftphase routine
c
        call zffti(nphi,wsave)

        do i=1,ngrid
        do j=1,nphi
        j0=j+nphi/2
        if( j0 .gt. nphi ) j0=j0-nphi
        grid_ctmp(j0)=grids_cmpl(i,j)
        enddo
ccc        call prin2('-grid_ctmp=*',grid_ctmp,2*nphi)
        call zfftb(nphi,grid_ctmp,wsave)
ccc        call prin2('+grid_ctmp=*',grid_ctmp,2*nphi)
        do j=1,nphi
        grids(i,j)=grid_ctmp(j)
        enddo
        enddo

        endif


        if( 1 .eq. 2 ) then
c
c       ... real version, needs shiftphase_real routine
c
        call dffti(nphi,wsave)

        do i=1,ngrid
        grid_rtmp(nphi)=0       
        j0=1
        j=j0+nphi/2
        grid_rtmp(j0)=dble(grids_cmpl(i,j))
        do j0=2,nphi/2
        j=j0+nphi/2
        grid_rtmp(2*j0-2)=+dble(grids_cmpl(i,j))
        grid_rtmp(2*j0-1)=+imag(grids_cmpl(i,j))
        enddo
ccc        call prin2('-grid_rtmp=*',grid_rtmp,nphi)
        call dfftb(nphi,grid_rtmp,wsave)
ccc        call prin2('+grid_rtmp=*',grid_rtmp,nphi)
        do j=1,nphi
        grids(i,j)=grid_rtmp(j)
        enddo
        enddo

        endif


        if( 2 .eq. 2 ) then
c
c       ... real version, needs shiftphase_real_alt routine
c
        call dffti(nphi,wsave)

        do i=1,ngrid
        grid_rtmp(nphi)=0       
        j0=1
        j=j0+nphi/2
        grid_rtmp(j0)=dble(grids_cmpl(i,j))
        do j0=2,nphi/2
        j=-j0+nphi/2+2
        grid_rtmp(2*j0-2)=+dble(grids_cmpl(i,j))
        grid_rtmp(2*j0-1)=-imag(grids_cmpl(i,j))
        enddo
ccc        call prin2('-grid_rtmp=*',grid_rtmp,nphi)
        call dfftb(nphi,grid_rtmp,wsave)
ccc        call prin2('+grid_rtmp=*',grid_rtmp,nphi)
        do j=1,nphi
        grids(i,j)=grid_rtmp(j)
        enddo
        enddo

        endif

        endif


        return
        end
c
c
c
c
c
        subroutine rotgrid_real_opt(nphi,phi,ntheta,theta,fmodes,
     $     ngrid,xgrid,ygrid,zgrid,grids)
        implicit real *8 (a-h,o-z)
c
c       Fast rotation of the user-defined grid (xgrid,ygrid,zgrid).
c       This function contructs a set of grids, obtained by rotating 
c       the original grid by nphi uniformly spaced angles phi.
c       
c       Assume symmetric spherical grid, nphi is even, ngrid = nphi*ntheta
c
c       The output is NGRID-by-NPHI real matrix, 
c       containing the rotated grids.
c
c       
c       Input parameters:
c       
c       nphi - the number of points in latitude discretization
c       phi - lattitude discretization angles
c       ntheta - the number of points in big circle discretization
c       theta - big circle discretization angles
c       fmodes - Fourier modes, NPHI-by-NTHETA complex matrix
c       ngrid - the number of points in the user-defined grid
c       xgrid - the x-coordinates of the grid
c       ygrid - the y-coordinates of the grid
c       zgrid - the z-coordinates of the grid
c       
c       Output parameters:
c       
c       grids - function values at rotated grids, NGRID-by-NPHI complex matrix
c
        real *8 phi(nphi),theta(ntheta)
        real *8 grids(ngrid,nphi)
        complex *16 fmodes(nphi,ntheta)
        real *8 xgrid(ngrid),ygrid(ngrid),zgrid(ngrid)
        real *8, allocatable :: zang(:)
        real *8, allocatable :: phase(:)
        real *8, allocatable :: phi_grid(:,:)
        real *8, allocatable :: theta_grid(:,:)

        complex *16, allocatable :: fmodes_tmp(:,:)
        complex *16, allocatable :: grid_ctmp(:)
        real *8, allocatable :: grid_rtmp(:)
        complex *16, allocatable :: wsave(:)

        complex *16, allocatable :: grids_cmpl(:,:)

        real *8, allocatable :: zang_tmp(:,:)
        complex *16, allocatable :: grids_tmp(:,:,:)

        eps=1d-12 * 1e0

c       storage for Fourier modes
ccc        lused = nphi*ntheta

c
c ...and the phases and angles to be interpolated are
c
        allocate(zang(ngrid))
        allocate(phase(ngrid))
c
ccc        lused = lused+2*ngrid
c
c
        do i=1,ngrid        
        call rotgrid3phitheta(xgrid(i),ygrid(i),zgrid(i),
     $     phase(i),zang(i))
        enddo
c
ccc        call prin2('zang=*',zang,ngrid)
ccc        call prin2('phase=*',phase,ngrid)


c Optionally, use nuFFT in R^2 to compare timings
c
        if_use_nufft2d = 0
        if( if_use_nufft2d .eq. 1 ) then
c
c Call nonuniform FFT in R^2 to get all grids simultaneously
c

        allocate( phi_grid(ngrid,nphi) )
        allocate( theta_grid(ngrid,nphi) )

        allocate( grids_cmpl(ngrid,nphi) )

        do j=1,nphi
        do i=1,ngrid
        phi_grid(i,j)=phi(j)+phase(i)
        theta_grid(i,j)=zang(i)
        enddo
        enddo

        iflag=+1
        call nufft2d2f90(ngrid*nphi,phi_grid,theta_grid,grids_cmpl,
     $     iflag,eps,nphi,ntheta,fmodes,ier)

        do j=1,nphi
        do i=1,ngrid
        grids(i,j)=grids_cmpl(i,j)
        enddo
        enddo

        return
        endif


        if_use_nufft1d = 1
        if( if_use_nufft1d .eq. 1 ) then
c
c Hybrid method
c Use nuFFT 1d in theta, and nuFFT/FFT 1d in phi
c

c
c Get interpolated values in theta
        iflag=+1
        allocate( fmodes_tmp(ntheta,nphi) )

        allocate( grids_cmpl(ngrid,nphi) )

ccc        lused = lused+2*ntheta*nphi + 2*ngrid*nphi

        if( 1 .eq. 2 ) then
c
c Reference nuFFT in R^1
c        do k=1,nphi
c        do i=1,ntheta
c        fmodes_tmp(i,k)=fmodes(k,i)
c        enddo
c        call nufft1d2f90(ngrid,zang,grids_cmpl(1,k),
c     $     iflag,eps,ntheta,fmodes_tmp(1,k),ier)
c        enddo

        endif
        
        if( 1 .eq. 2 ) then
c
c Vectorized nuFFT in R^1
c        do k=1,nphi
c        do i=1,ntheta
c        fmodes_tmp(i,k)=fmodes(k,i)
c        enddo
c        enddo
c        call nufft1d2vf90(nphi,ngrid,zang,grids_cmpl,
c     $     iflag,eps,ntheta,fmodes_tmp,ier)

c        do i=1,nphi
c        call prinf('i=*',i,1)
c        call prin2('grids_cmpl=*',grids_cmpl(1,i),2*ngrid)
c        enddo
c        stop
        
        endif

        if( 1 .eq. 2 ) then
c
c Vectorized nuFFT in R^1, with real symmetries
        do k=1,nphi/2+1
        do i=1,ntheta
        fmodes_tmp(i,k)=fmodes(k,i)
        enddo
        enddo

c        do i=1,nphi/2
c        call prinf('i=*',i,1)
c        call prin2('fmodes_tmp=*',fmodes_tmp(1,i),2*ntheta)
c        enddo
c        stop

        call nufft1d2vf90(nphi/2+1,ngrid,zang,grids_cmpl,
     $     iflag,eps,ntheta,fmodes_tmp,ier)

c        write(*,*) ngrid, ntheta, nphi
c        do i=1,nphi/2+1
c        call prinf('i=*',i,1)
c        call prin2('grids_cmpl+=*',grids_cmpl(1,i),2*ngrid)
c        enddo
c        stop

        endif

        if( 2 .eq. 2 ) then
c
c Assume symmetric grid, nphi is even        
c Vectorized nuFFT in R^1, with real symmetries
        do k=1,nphi/2+1
        do i=1,ntheta
        fmodes_tmp(i,k)=fmodes(k,i)
        enddo
        enddo

        allocate(zang_tmp(nphi/2+1,ntheta/2))
        allocate(grids_tmp(nphi/2+1,ntheta/2,nphi/2+1))

ccc        lused = lused+(nphi/2+1)*ntheta/2
ccc        lused = lused+2*(nphi/2+1)*(ntheta/2)*(nphi/2+1)

        kk=1
        do i=1,ntheta/2
        do k=1,nphi/2+1
        zang_tmp(k,i)=zang(kk)
        kk=kk+1
        enddo
        kk=1+i*nphi
        enddo

c        do i=1,nphi/2
c        call prinf('i=*',i,1)
c        call prin2('fmodes_tmp=*',fmodes_tmp(1,i),2*ntheta)
c        enddo
c        stop

        ngrid_tmp = (nphi/2+1)*ntheta/2
        call nufft1d2vf90(nphi/2+1,ngrid_tmp,zang_tmp,grids_tmp,
     $     iflag,eps,ntheta,fmodes_tmp,ier)


c        write(*,*) ngrid_tmp, ntheta, nphi
c        do i=1,nphi/2+1
c        call prinf('i=*',i,1)
c        call prin2('grids_tmp=*',grids_tmp(1,1,i),2*ngrid_tmp)
c        enddo

        do j=1,nphi/2+1
        kk=1
        do i=1,ntheta/2
        do k=1,nphi/2+1
        grids_cmpl(kk,j)=grids_tmp(k,i,j)
        kk=kk+1
        enddo
        do k=nphi/2+2,nphi
        grids_cmpl(kk,j)=grids_tmp(nphi+2-k,i,j)
        kk=kk+1
        enddo
        enddo
        enddo

c        write(*,*) ngrid, ntheta, nphi
c        do i=1,nphi/2+1
c        call prinf('i=*',i,1)
c        call prin2('grids_cmpl=*',grids_cmpl(1,i),2*ngrid)
c        enddo
c        stop

        endif



        if( 1 .eq. 2 ) then
c
c initialize the fourier grid for shiftphase and shiftphase_real routine
c
        do i=nphi/2+2,nphi
cc        write(*,*) i,nphi-i+2
        i0=nphi-i+2
        do j=1,ngrid
        grids_cmpl(j,i)=dconjg(grids_cmpl(j,i0))
        enddo
        enddo

        endif

c
c Shift the grids on latitudes (apply phase in Fourier domain)
ccc        call shiftphase(grids_cmpl,ngrid,nphi,phase)
ccc        call shiftphase_real(grids_cmpl,ngrid,nphi,phase)
        call shiftphase_real_alt(grids_cmpl,ngrid,nphi,phase)


c
c FFT data on latitudes, apply inverse fftshift, even nphi only
c octave:9> ifftshift( 1:5 )
c ans =
c
c    3   4   5   1   2
c
c octave:10> ifftshift( 1:6 )
c ans =
c
c    4   5   6   1   2   3

        allocate( grid_ctmp(nphi) )
        allocate( grid_rtmp(nphi) )
        allocate( wsave(4*nphi+15) )


ccc        call prinf('lused=*',lused,1)
ccc        call prinf('lused(k)=*',lused/1000,1)
ccc        call prinf('lused(m)=*',lused/1000000,1)



        if ( 1 .eq. 2 ) then
c
c       ... complex version, needs shiftphase routine
c
        call zffti(nphi,wsave)

        do i=1,ngrid
        do j=1,nphi
        j0=j+nphi/2
        if( j0 .gt. nphi ) j0=j0-nphi
        grid_ctmp(j0)=grids_cmpl(i,j)
        enddo
ccc        call prin2('-grid_ctmp=*',grid_ctmp,2*nphi)
        call zfftb(nphi,grid_ctmp,wsave)
ccc        call prin2('+grid_ctmp=*',grid_ctmp,2*nphi)
        do j=1,nphi
        grids(i,j)=grid_ctmp(j)
        enddo
        enddo

        endif


        if( 1 .eq. 2 ) then
c
c       ... real version, needs shiftphase_real routine
c
        call dffti(nphi,wsave)

        do i=1,ngrid
        grid_rtmp(nphi)=0       
        j0=1
        j=j0+nphi/2
        grid_rtmp(j0)=dble(grids_cmpl(i,j))
        do j0=2,nphi/2
        j=j0+nphi/2
        grid_rtmp(2*j0-2)=+dble(grids_cmpl(i,j))
        grid_rtmp(2*j0-1)=+imag(grids_cmpl(i,j))
        enddo
ccc        call prin2('-grid_rtmp=*',grid_rtmp,nphi)
        call dfftb(nphi,grid_rtmp,wsave)
ccc        call prin2('+grid_rtmp=*',grid_rtmp,nphi)
        do j=1,nphi
        grids(i,j)=grid_rtmp(j)
        enddo
        enddo

        endif


        if( 2 .eq. 2 ) then
c
c       ... real version, needs shiftphase_real_alt routine
c
        call dffti(nphi,wsave)

        do i=1,ngrid
        grid_rtmp(nphi)=0       
        j0=1
        j=j0+nphi/2
        grid_rtmp(j0)=dble(grids_cmpl(i,j))
        do j0=2,nphi/2
        j=-j0+nphi/2+2
        grid_rtmp(2*j0-2)=+dble(grids_cmpl(i,j))
        grid_rtmp(2*j0-1)=-imag(grids_cmpl(i,j))
        enddo
ccc        call prin2('-grid_rtmp=*',grid_rtmp,nphi)
        call dfftb(nphi,grid_rtmp,wsave)
ccc        call prin2('+grid_rtmp=*',grid_rtmp,nphi)
        do j=1,nphi
        grids(i,j)=grid_rtmp(j)
        enddo
        enddo

        endif

        endif


        return
        end
c
c
c
c
c
        subroutine rotgrid_real_ref(nphi,phi,ntheta,theta,fmodes,
     $     ngrid,xgrid,ygrid,zgrid,grids)
        implicit real *8 (a-h,o-z)
c
c       Fast rotation of the user-defined grid (xgrid,ygrid,zgrid).
c       This function contructs a set of grids, obtained by rotating 
c       the original grid by nphi uniformly spaced angles phi.
c       
c       Assume symmetric spherical grid, nphi is even, ngrid = nphi*ntheta
c
c       The output is NGRID-by-NPHI real matrix, 
c       containing the rotated grids.
c
c       
c       Input parameters:
c       
c       nphi - the number of points in latitude discretization
c       phi - lattitude discretization angles
c       ntheta - the number of points in big circle discretization
c       theta - big circle discretization angles
c       fmodes - Fourier modes, NPHI-by-NTHETA complex matrix
c       ngrid - the number of points in the user-defined grid
c       xgrid - the x-coordinates of the grid
c       ygrid - the y-coordinates of the grid
c       zgrid - the z-coordinates of the grid
c       
c       Output parameters:
c       
c       grids - function values at rotated grids, NGRID-by-NPHI complex matrix
c
        real *8 phi(nphi),theta(ntheta)
        real *8 grids(ngrid,nphi)
        complex *16 fmodes(nphi,ntheta)
        real *8 xgrid(ngrid),ygrid(ngrid),zgrid(ngrid)
        real *8, allocatable :: zang(:)
        real *8, allocatable :: phase(:)
        real *8, allocatable :: phi_grid(:,:)
        real *8, allocatable :: theta_grid(:,:)

        complex *16, allocatable :: fmodes_tmp(:,:)
        complex *16, allocatable :: grid_ctmp(:)
        real *8, allocatable :: grid_rtmp(:)
        complex *16, allocatable :: wsave(:)

        complex *16, allocatable :: grids_cmpl(:,:)

        real *8, allocatable :: zang_tmp(:,:)
        complex *16, allocatable :: grids_tmp(:,:,:)

        eps=1d-12 * 1e0

c
c ...and the phases and angles to be interpolated are
c
        allocate(zang(ngrid))
        allocate(phase(ngrid))
c
        do i=1,ngrid        
        call rotgrid3phitheta(xgrid(i),ygrid(i),zgrid(i),
     $     phase(i),zang(i))
        enddo
c
ccc        call prin2('zang=*',zang,ngrid)
ccc        call prin2('phase=*',phase,ngrid)


c
c Optionally, use nuFFT in R^2 to compare timings
c
        if_use_nufft2d = 1
        if( if_use_nufft2d .eq. 1 ) then
c
c Call nonuniform FFT in R^2 to get all grids simultaneously
c

        allocate( phi_grid(ngrid,nphi) )
        allocate( theta_grid(ngrid,nphi) )

        allocate( grids_cmpl(ngrid,nphi) )

        do j=1,nphi
        do i=1,ngrid
        phi_grid(i,j)=phi(j)+phase(i)
        theta_grid(i,j)=zang(i)
        enddo
        enddo

        iflag=+1
        call nufft2d2f90(ngrid*nphi,phi_grid,theta_grid,grids_cmpl,
     $     iflag,eps,nphi,ntheta,fmodes,ier)

        do j=1,nphi
        do i=1,ngrid
        grids(i,j)=grids_cmpl(i,j)
        enddo
        enddo

        return
        endif

        return
        end
c
c
c
c
c
        subroutine rotgrid_real_kbw(nphi,phi,ntheta,theta,fmodes,
     $     ngrid,xgrid,ygrid,zgrid,grids,w)
        implicit real *8 (a-h,o-z)
c
c       Fast rotation of the user-defined grid (xgrid,ygrid,zgrid).
c       This function contructs a set of grids, obtained by rotating 
c       the original grid by nphi uniformly spaced angles phi.
c       
c       Assume symmetric spherical grid, nphi is even, ngrid = nphi*ntheta
c
c       The output is NGRID-by-NPHI real matrix, 
c       containing the rotated grids.
c
c       
c       Input parameters:
c       
c       nphi - the number of points in latitude discretization
c       phi - lattitude discretization angles
c       ntheta - the number of points in big circle discretization
c       theta - big circle discretization angles
c       fmodes - Fourier modes, NPHI-by-NTHETA complex matrix
c       ngrid - the number of points in the user-defined grid
c       xgrid - the x-coordinates of the grid
c       ygrid - the y-coordinates of the grid
c       zgrid - the z-coordinates of the grid
c       
c       Output parameters:
c       
c       grids - function values at rotated grids, NGRID-by-NPHI complex matrix
c
        real *8 phi(nphi),theta(ntheta)
        real *8 grids(ngrid,nphi)
        complex *16 fmodes(nphi,ntheta)
        real *8 xgrid(ngrid),ygrid(ngrid),zgrid(ngrid)
        real *8, allocatable :: zang(:)
        real *8, allocatable :: phase(:)
        real *8, allocatable :: phi_grid(:,:)
        real *8, allocatable :: theta_grid(:,:)

        complex *16, allocatable :: fmodes_tmp(:,:)
        complex *16, allocatable :: grid_ctmp(:)
        real *8, allocatable :: grid_rtmp(:)
        complex *16, allocatable :: wsave(:)

        complex *16, allocatable :: grids_cmpl(:,:)

        real *8, allocatable :: zang_tmp(:,:)
        complex *16, allocatable :: grids_tmp(:,:,:)

        real *8 w(*)

        eps=1d-12 * 1e0

c
c ...and the phases and angles to be interpolated are
c
        allocate(zang(ngrid))
        allocate(phase(ngrid))
c
        do i=1,ngrid        
        call rotgrid3phitheta(xgrid(i),ygrid(i),zgrid(i),
     $     phase(i),zang(i))
        enddo
c
ccc        call prin2('zang=*',zang,ngrid)
ccc        call prin2('phase=*',phase,ngrid)


c
c Optionally, use nuFFT in R^2 to compare timings
c
        if_use_nufft2d = 1
        if( if_use_nufft2d .eq. 1 ) then
c
c Call nonuniform FFT in R^2 to get all grids simultaneously
c

        allocate( phi_grid(ngrid,nphi) )
        allocate( theta_grid(ngrid,nphi) )

        allocate( grids_cmpl(ngrid,nphi) )

        do j=1,nphi
        do i=1,ngrid
        phi_grid(i,j)=phi(j)+phase(i)
        theta_grid(i,j)=zang(i)
        enddo
        enddo

        iflag=+1
        call nufft2d2f90kbw(ngrid*nphi,phi_grid,theta_grid,grids_cmpl,
     $     iflag,eps,nphi,ntheta,fmodes,ier,w)

        do j=1,nphi
        do i=1,ngrid
        grids(i,j)=grids_cmpl(i,j)
        enddo
        enddo

        return
        endif

        return
        end
c
c
c
c
c
        subroutine rotgrid_real_kb0(nphi,phi,ntheta,theta,fmodes,
     $     ngrid,xgrid,ygrid,zgrid,grids,lused)
        implicit real *8 (a-h,o-z)
c
c       Fast rotation of the user-defined grid (xgrid,ygrid,zgrid).
c       This function contructs a set of grids, obtained by rotating 
c       the original grid by nphi uniformly spaced angles phi.
c       
c       Assume symmetric spherical grid, nphi is even, ngrid = nphi*ntheta
c
c       The output is NGRID-by-NPHI real matrix, 
c       containing the rotated grids.
c
c       
c       Input parameters:
c       
c       nphi - the number of points in latitude discretization
c       phi - lattitude discretization angles
c       ntheta - the number of points in big circle discretization
c       theta - big circle discretization angles
c       fmodes - Fourier modes, NPHI-by-NTHETA complex matrix
c       ngrid - the number of points in the user-defined grid
c       xgrid - the x-coordinates of the grid
c       ygrid - the y-coordinates of the grid
c       zgrid - the z-coordinates of the grid
c       
c       Output parameters:
c       
c       grids - function values at rotated grids, NGRID-by-NPHI complex matrix
c
        real *8 phi(nphi),theta(ntheta)
        real *8 grids(ngrid,nphi)
        complex *16 fmodes(nphi,ntheta)
        real *8 xgrid(ngrid),ygrid(ngrid),zgrid(ngrid)
        real *8, allocatable :: zang(:)
        real *8, allocatable :: phase(:)
        real *8, allocatable :: phi_grid(:,:)
        real *8, allocatable :: theta_grid(:,:)

        complex *16, allocatable :: fmodes_tmp(:,:)
        complex *16, allocatable :: grid_ctmp(:)
        real *8, allocatable :: grid_rtmp(:)
        complex *16, allocatable :: wsave(:)

        complex *16, allocatable :: grids_cmpl(:,:)

        real *8, allocatable :: zang_tmp(:,:)
        complex *16, allocatable :: grids_tmp(:,:,:)

        eps=1d-12 * 1e0

c
c ...and the phases and angles to be interpolated are
c
        allocate(zang(ngrid))
        allocate(phase(ngrid))
c
        do i=1,ngrid        
        call rotgrid3phitheta(xgrid(i),ygrid(i),zgrid(i),
     $     phase(i),zang(i))
        enddo
c
ccc        call prin2('zang=*',zang,ngrid)
ccc        call prin2('phase=*',phase,ngrid)


c
c Optionally, use nuFFT in R^2 to compare timings
c
        if_use_nufft2d = 1
        if( if_use_nufft2d .eq. 1 ) then
c
c Call nonuniform FFT in R^2 to get all grids simultaneously
c

        iflag=+1
        call nufft2dkb0(ngrid*nphi,eps,nphi,ntheta,lused)

        return
        endif

        return
        end
c
c
c
c
c
        subroutine rotgrid_real_kbi(nphi,phi,ntheta,theta,fmodes,
     $     ngrid,xgrid,ygrid,zgrid,grids,w,lw,lused)
        implicit real *8 (a-h,o-z)
c
c       Fast rotation of the user-defined grid (xgrid,ygrid,zgrid).
c       This function contructs a set of grids, obtained by rotating 
c       the original grid by nphi uniformly spaced angles phi.
c       
c       Assume symmetric spherical grid, nphi is even, ngrid = nphi*ntheta
c
c       The output is NGRID-by-NPHI real matrix, 
c       containing the rotated grids.
c
c       
c       Input parameters:
c       
c       nphi - the number of points in latitude discretization
c       phi - lattitude discretization angles
c       ntheta - the number of points in big circle discretization
c       theta - big circle discretization angles
c       fmodes - Fourier modes, NPHI-by-NTHETA complex matrix
c       ngrid - the number of points in the user-defined grid
c       xgrid - the x-coordinates of the grid
c       ygrid - the y-coordinates of the grid
c       zgrid - the z-coordinates of the grid
c       
c       Output parameters:
c       
c       grids - function values at rotated grids, NGRID-by-NPHI complex matrix
c
        real *8 phi(nphi),theta(ntheta)
        real *8 grids(ngrid,nphi)
        complex *16 fmodes(nphi,ntheta)
        real *8 xgrid(ngrid),ygrid(ngrid),zgrid(ngrid)
        real *8, allocatable :: zang(:)
        real *8, allocatable :: phase(:)
        real *8, allocatable :: phi_grid(:,:)
        real *8, allocatable :: theta_grid(:,:)

        complex *16, allocatable :: fmodes_tmp(:,:)
        complex *16, allocatable :: grid_ctmp(:)
        real *8, allocatable :: grid_rtmp(:)
        complex *16, allocatable :: wsave(:)

        complex *16, allocatable :: grids_cmpl(:,:)

        real *8, allocatable :: zang_tmp(:,:)
        complex *16, allocatable :: grids_tmp(:,:,:)

        real *8 w(*)

        eps=1d-12 * 1e0

c
c ...and the phases and angles to be interpolated are
c
        allocate(zang(ngrid))
        allocate(phase(ngrid))
c
        do i=1,ngrid        
        call rotgrid3phitheta(xgrid(i),ygrid(i),zgrid(i),
     $     phase(i),zang(i))
        enddo
c
ccc        call prin2('zang=*',zang,ngrid)
ccc        call prin2('phase=*',phase,ngrid)


c
c Optionally, use nuFFT in R^2 to compare timings
c
        if_use_nufft2d = 1
        if( if_use_nufft2d .eq. 1 ) then
c
c Call nonuniform FFT in R^2 to get all grids simultaneously
c

        allocate( phi_grid(ngrid,nphi) )
        allocate( theta_grid(ngrid,nphi) )

        allocate( grids_cmpl(ngrid,nphi) )

        do j=1,nphi
        do i=1,ngrid
        phi_grid(i,j)=phi(j)+phase(i)
        theta_grid(i,j)=zang(i)
        enddo
        enddo

        iflag=+1
        call nufft2dkb0(ngrid*nphi,eps,nphi,ntheta,lused)
        call nufft2dkbi(ngrid*nphi,phi_grid,theta_grid,
     $     eps,nphi,ntheta,ier,w,lw,lused)

        return
        endif

        return
        end
c
c
c
c
c
