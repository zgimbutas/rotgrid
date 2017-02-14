addpath '../matlab';
addpath '../nufft1d';
addpath '../nufft2d';


nterms = 31
beta = pi/3

nrot = 2*nterms+2
nrot = fftnext235(nrot)

mpole = zeros(nterms+1,2*nterms+1)+1i*zeros(nterms+1,2*nterms+1);

jc = nterms+1;
for j=1:nterms+1
  mpole(j,jc+(-j+1:j-1)) = 1; 
end

j=nterms+1;
for i=0:nterms
mpole(i+1,j+(-i:i)) = rand(1,2*i+1) + 1i*rand(1,2*i+1);
end

mpole;


'Wigner rotation matrices + FFT'
ldc = nterms;
tic
rotmat = real(rotmat_proj_init(ldc,beta));
toc
tic
mpout=rot1lat_fsr_cmpl(nterms,mpole,beta,nrot,rotmat,ldc);
toc
mpout=reshape(mpout,nterms+1,2*nterms+1,nrot);


'Pseudo-spectral projection scheme'
tic
mptmp=rot1lat_proj_cmpl(nterms,mpole,beta,nrot);
toc
mptmp=reshape(mptmp,nterms+1,2*nterms+1,nrot);


error = norm(reshape(mpout-mptmp,(nterms+1)*(2*nterms+1),nrot),2)


