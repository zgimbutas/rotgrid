
nterms = 48

nrot = 2*nterms+2
nrot = fftnext235(nrot)

mpole = zeros(nterms+1,2*nterms+1)+1i*zeros(nterms+1,2*nterms+1);

jc = nterms+1;
for j=1:nterms+1
  mpole(j,jc+(-j+1:j-1)) = 1; 
end

mpole;

itype = 2;

if( itype == 2 ),

ldc = nterms;
nbeta = nterms+1;

beta = grule(nbeta);
rotmat = real(rotmat_proj_init(ldc,beta(1)));
[n,m] = size(rotmat)
rotmat = zeros(n,nbeta);
tic
for i = 1:nbeta
  rotmat(:,i) = real(rotmat_proj_init(ldc,beta(i)));
end
toc

sKwords = ( n * nbeta ) /1e3
sMwords = ( n * nbeta ) /1e6
sGwords = ( n * nbeta ) /1e9

mpout = zeros(nterms+1,2*nterms+1,nrot,nbeta);
tic
for i = 1:nbeta
mptmp=rot1lat_fsr_cmpl(nterms,mpole,beta(i),nrot,rotmat(:,i),ldc);
mpout(:,:,:,i)=reshape(mptmp,nterms+1,2*nterms+1,nrot);
end
toc
tic
mpout=reshape(mpout,nterms+1,2*nterms+1,nrot*nbeta);
toc

sKwords = ((nterms+1) * (2*nterms+1) * nrot*nbeta) /1e3
sMwords = ((nterms+1) * (2*nterms+1) * nrot*nbeta) /1e6
sGwords = ((nterms+1) * (2*nterms+1) * nrot*nbeta) /1e9

end

