%  Testing routines for fast gridding algorithm
%
%  Fast rotation of the user-defined grid (xrot,yrot,zrot).
%  This function constructs a set of grids, obtained by rotating 
%  the original grid by nphi uniformly spaced angles phi.


nterms = 12*4
mpole = zeros(nterms+1,2*nterms+1)+1i*zeros(nterms+1,2*nterms+1);

%j=nterms+1;
%for i=0:nterms
%mpole(i+1,j+(-i:i)) = 1/sqrt(2*i+1); 
%end

j=nterms+1;
for i=0:nterms
mpole(i+1,j+(-i:i)) = rand(1,2*i+1) + 1i*rand(1,2*i+1);
end

ntheta = nterms+1
nphi = 2*nterms+2
nphi = fftnext235(nphi)
ngrid = ntheta*nphi

[ctheta,whts,ynms,wsave]=sphtrans_cmpl_cheb_init(nterms,nphi,ntheta);

fgrid = sphtrans_cmpl(nterms,mpole,nphi,ntheta,ctheta,ynms,wsave);

% double grid
a = fgrid(1:nphi/2,:);
b = fgrid(nphi/2+1:nphi,:);
fgrid = [[a; b] [fliplr(b); fliplr(a)]];


beta =  acos(ctheta);
nbeta = ntheta;

nrot = nphi;

[rotmat]=rotgrid_fsr_cmpl_init(nterms,nbeta,beta);

tic
grids=rotgrid_dsr_cmpl(nterms,mpole,nphi,ntheta,nbeta,beta,nrot,...
           rotmat,ctheta,ynms,wsave);
toc

grids=reshape(grids,nphi,ntheta,nrot,nbeta);

%
%  Test the rotated grids
%
k=2;

for j=1:nrot
mpout = rotviaproj_cmpl(nterms,mpole,beta(k),+2*pi*(j-1)/nrot);
ftest = sphtrans_cmpl(nterms,mpout,nphi,ntheta,ctheta,ynms,wsave);
errors(j)=norm(grids(:,:,j,k)-ftest,2);
end

norm(errors/nrot)
