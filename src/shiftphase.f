
        subroutine shiftphase(amatr,ngrid,nphi,phase)
        implicit real *8 (a-h,o-z)
c
c       This routine evaluates products
c       amatr(i,:) = amatr(i,:) .* exp(+1i*((-nphi/2:(nphi-1)/2))*phase(i));
c
        complex *16 amatr(ngrid,nphi),ima,ztmp
        real *8 phase(ngrid)
        complex *16, allocatable :: ephi0(:)
        data ima/(0.0d0,1.0d0)/
c
c
        if( 2 .eq. 2 ) then
c
        allocate( ephi0(-nphi/2:+nphi/2) )
c
        do i=1,ngrid
c
        ztmp=exp(+ima*phase(i))
        ephi0(0)=1
        ephi0(1)=ztmp
        ephi0(-1)=dconjg(ztmp)
        do j=2,nphi/2
        ephi0(j)=ephi0(j-1)*ztmp
        ephi0(-j)=dconjg(ephi0(j))
        enddo
c
        do j=-nphi/2,(nphi-1)/2
        amatr(i,j+nphi/2+1) = amatr(i,j+nphi/2+1)*ephi0(j)
        enddo
c
        enddo
c
        endif
c
c
c
        if( 1 .eq. 2 ) then

        do i=1,ngrid
c
        do j=-nphi/2,(nphi-1)/2
        amatr(i,j+nphi/2+1) = amatr(i,j+nphi/2+1)*exp(+ima*j*phase(i))
        enddo
c
        enddo
c
        endif
c
c
        return
        end
c
c
c
c
        subroutine shiftphase_real(amatr,ngrid,nphi,phase)
        implicit real *8 (a-h,o-z)
c
c       This routine evaluates products
c       amatr(i,:) = amatr(i,:) .* exp(+1i*((-nphi/2:(nphi-1)/2))*phase(i));
c
c       Shift positive frequencies only, for range nphi/2+1:nphi
c
        complex *16 amatr(ngrid,nphi),ima,ztmp
        real *8 phase(ngrid)
        complex *16, allocatable :: ephi0(:)
        data ima/(0.0d0,1.0d0)/
c
c
        if( 2 .eq. 2 ) then
c
        allocate( ephi0(0:+nphi/2) )
c
        do i=1,ngrid
c
        ztmp=exp(+ima*phase(i))
        ephi0(0)=1
        ephi0(1)=ztmp
        do j=2,nphi/2
        ephi0(j)=ephi0(j-1)*ztmp
        enddo
c
        do j=0,(nphi-1)/2
        amatr(i,j+nphi/2+1) = amatr(i,j+nphi/2+1)*ephi0(j)
        enddo
c
        enddo
c
        endif
c
c
        return
        end
c
c
c
c
c
        subroutine shiftphase_real_alt(amatr,ngrid,nphi,phase)
        implicit real *8 (a-h,o-z)
c
c       This routine evaluates products
c       amatr(i,:) = amatr(i,:) .* exp(+1i*((-nphi/2:(nphi-1)/2))*phase(i));
c
c       Shift negative frequencies only, for range 1:nphi/2+1
c
        complex *16 amatr(ngrid,nphi),ima,ztmp
        real *8 phase(ngrid)
        complex *16, allocatable :: ephi0(:)
        data ima/(0.0d0,1.0d0)/
c
c
        if( 2 .eq. 2 ) then
c
        allocate( ephi0(0:+nphi/2) )
c
        do i=1,ngrid
c
        ztmp=exp(-ima*phase(i))
        ephi0(0)=1
        ephi0(1)=ztmp
        do j=2,nphi/2
        ephi0(j)=ephi0(j-1)*ztmp
        enddo
c
        do j=0,(nphi-1)/2
        amatr(i,-j+nphi/2+1) = amatr(i,-j+nphi/2+1)*ephi0(j)
        enddo
c
        enddo
c
        endif
c
c
        return
        end
c
c
c
c
c
