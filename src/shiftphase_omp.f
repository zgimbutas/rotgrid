
        subroutine shiftphase(cii,nii,nphi,pii)
        implicit real *8 (a-h,o-z)
c
c       This routine evaluates products
c       cii(i,:) = cii(i,:) .* exp(+1i*((-nphi/2:(nphi-1)/2))*pii(i));
c
        complex *16 cii(nii,nphi),ima,ztmp
        real *8 pii(nii)
        complex *16, allocatable :: ephi0(:)
        data ima/(0.0d0,1.0d0)/
c
c
        if( 2 .eq. 2 ) then
c
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ztmp,ephi0,i,j)
        do i=1,nii
c
        allocate( ephi0(-nphi/2:+nphi/2) )
c
        ztmp=exp(+ima*pii(i))
        ephi0(0)=1
        ephi0(1)=ztmp
        ephi0(-1)=dconjg(ztmp)
        do j=2,nphi/2
        ephi0(j)=ephi0(j-1)*ztmp
        ephi0(-j)=dconjg(ephi0(j))
        enddo
c
        do j=-nphi/2,(nphi-1)/2
        cii(i,j+nphi/2+1) = cii(i,j+nphi/2+1)*ephi0(j)
        enddo
c
        deallocate(ephi0)
c
        enddo
C$OMP END PARALLEL DO
c
        endif
c
c
c
        if( 1 .eq. 2 ) then

        do i=1,nii
c
        do j=-nphi/2,(nphi-1)/2
        cii(i,j+nphi/2+1) = cii(i,j+nphi/2+1)*exp(+ima*j*pii(i))
        enddo
c
        enddo
c
        endif
c
c
        return
        end
