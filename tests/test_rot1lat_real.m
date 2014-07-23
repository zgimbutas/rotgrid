
nterms = 31
beta = pi/3

nrot = 2*nterms+2
nrot = fftnext235(nrot)

mpole = zeros(nterms+1,nterms+1)+1i*zeros(nterms+1,nterms+1);

jc = nterms+1;
for j=1:nterms+1
  mpole(j,jc+(1:j-1)) = 1; 
end

j=nterms+1;
for i=0:nterms
mpole(i+1,j+(1:i)) = rand(1,i) + 1i*rand(1,i);
end

mpole;


'Wigner rotation matrices + FFT'
ldc = nterms;
tic
rotmat = real(rotmat_proj_init(ldc,beta));
toc
tic
mpout=rot1lat_fsr_real(nterms,mpole,beta,nrot,rotmat,ldc);
toc
mpout=reshape(mpout,nterms+1,nterms+1,nrot);

