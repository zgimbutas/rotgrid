c
c       Vectorized nuFFT1d routines + OpenMP
c

        subroutine nufft1d1vf90(nvec,nj,xj,cj,iflag,eps,ms,fk,ier)
        implicit real *8 (a-h,o-z)
        real*8 xj(nj)
        complex*16 cj(nj,nvec), fk(-ms/2:(ms-1)/2,nvec)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k)
        do k=1,nvec
        call nufft1d1f90(nj,xj,cj(1,k),iflag,eps,ms,fk(-ms/2,k),ier)
        enddo
C$OMP END PARALLEL DO

        return
        end


        subroutine nufft1d2vf90(nvec,nj,xj,cj,iflag,eps,ms,fk,ier)
        implicit real *8 (a-h,o-z)
        real*8 xj(nj)
        complex*16 cj(nj,nvec), fk(-ms/2:(ms-1)/2,nvec)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k)
        do k=1,nvec
        call nufft1d2f90(nj,xj,cj(1,k),iflag,eps,ms,fk(-ms/2,k),ier)
        enddo
C$OMP END PARALLEL DO

        return
        end


        subroutine nufft1d3vf90(nvec,nj,xj,cj,iflag,eps,nk,sk,fk,ier)
        implicit real *8 (a-h,o-z)
        real*8 xj(nj), sk(nk)
        complex*16 cj(nj,nvec), fk(nk,nvec)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k)
        do k=1,nvec
        call nufft1d3f90(nj,xj,cj(1,k),iflag,eps,nk,sk,fk(1,k),ier)
        enddo
C$OMP END PARALLEL DO

        return
        end
