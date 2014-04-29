%  Testing routines for spherical transform algorithm
%

nterms = 32
mpole = zeros(nterms+1,nterms+1)+1i*zeros(nterms+1,nterms+1);

j=1;
for i=0:nterms
mpole(i+1,j+(1:i)) = rand(1,i) + 1i*rand(1,i);
end


ntheta = nterms+1
nphi = 2*nterms+2
nphi = fftnext235(nphi)

[ctheta,whts,ynms,wsave]=sphtrans_real_lege_init(nterms,nphi,ntheta);

fgrid = sphtrans_real(nterms,mpole,nphi,ntheta,ctheta,ynms,wsave);
mpout = sphtrans_fwd_real(nterms,nphi,ntheta,fgrid,ctheta,whts,ynms,wsave);

error=norm(mpole-mpout,2)
