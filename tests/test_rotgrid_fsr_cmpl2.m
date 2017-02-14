%  Testing routines for fast gridding algorithm
%
%  Fast rotation of the user-defined grid (xrot,yrot,zrot).
%  This function constructs a set of grids, obtained by rotating 
%  the original grid by nphi uniformly spaced angles phi.

addpath '../matlab';
addpath '../nufft1d';
addpath '../nufft2d';


p = 25;
nterms = p+1;

ntheta=nterms;
nphi=2*nterms+2; 
nphi=fftnext235(nphi);

%
% Get the values of the user-defined function 
% on the Legendre grid (theta in [0..pi], see note above)
%
[phi,theta,xs,ys,zs]=init_grid_lege_single(nphi,ntheta);
beta=0;
[xgrid,ygrid,zgrid]=xyz_grid(beta,nphi,xs,ys,ntheta,zs,theta);
fgrid = funuser(xgrid,ygrid,zgrid);

%
% Evaluate the coefficients of spherical harmonic expansion
%
[ctheta,whts,ynms,wsave]=sphtrans_cmpl_lege_init(nterms,nphi,ntheta);
mpole = sphtrans_fwd_cmpl(nterms,nphi,ntheta,fgrid,ctheta,whts,ynms,wsave);


%
% Evaluate function values at all rotate grids
%
beta =  acos(ctheta);
nbeta = ntheta;

nrot = nphi;

[rotmat]=rotgrid_fsr_cmpl_init(nterms,nbeta,beta);

tic
grids=rotgrid_fsr_cmpl(nterms,mpole,nphi,ntheta,nbeta,beta,nrot,...
           rotmat,ctheta,ynms,wsave);
toc

grids=reshape(grids,nphi,ntheta,nrot,nbeta);


%
%  Test the rotated grids via rotated spherical harmonic expansions
%
k=2;

errors=zeros(1,nphi);
for j=1:nrot
mpout = rotviaproj_cmpl(nterms,mpole,beta(k),+2*pi*(j-1)/nrot);
ftest = sphtrans_cmpl(nterms,mpout,nphi,ntheta,ctheta,ynms,wsave);
errors(j)=norm(grids(:,:,j,k)-ftest,2);
end

%%%errors
max(abs(errors))


% 
% Finally, check the errors for all rotations directly
%
k=2;

errors=zeros(1,nrot);
for j=1:nrot
alpha=2*pi*(j-1)/nrot;
[xrota,yrota,zrota]=xyz_grid_a(beta(k),nphi,xs,ys,ntheta,zs,theta,alpha);
ftest = funuser(xrota,yrota,zrota);
errors(j)=norm(grids(:,:,j,k)-ftest,2);
end

%%%errors
max(abs(errors))
