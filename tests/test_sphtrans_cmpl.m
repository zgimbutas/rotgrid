%  Testing routines for spherical transform algorithm
%

addpath '../matlab';
addpath '../nufft1d';
addpath '../nufft2d';


nterms = 32
mpole = zeros(nterms+1,2*nterms+1)+1i*zeros(nterms+1,2*nterms+1);

j=nterms+1;
for i=0:nterms
mpole(i+1,j+(-i:i)) = rand(1,2*i+1) + 1i*rand(1,2*i+1);
end

ntheta = nterms+1
nphi = 2*nterms+2
nphi = fftnext235(nphi)

[ctheta,whts,ynms,wsave]=sphtrans_cmpl_lege_init(nterms,nphi,ntheta);

fgrid = sphtrans_cmpl(nterms,mpole,nphi,ntheta,ctheta,ynms,wsave);
mpout = sphtrans_fwd_cmpl(nterms,nphi,ntheta,fgrid,ctheta,whts,ynms,wsave);

error=norm(mpole-mpout,2)
