%  Testing routines for fast gridding algorithm
%
%  Fast rotation of the user-defined grid (xrot,yrot,zrot).
%  This function constructs a set of grids, obtained by rotating 
%  the original grid by nphi uniformly spaced angles phi.

addpath '../src';
addpath '../nufft1d';
addpath '../nufft2d';


%
% Construct the initial spherical grid, 
% it will be used to get Fourier coefficients for a sphere.
%
% NOTE: we are using great circles theta = [-pi..pi] 
% (and not the meridians theta = [0..pi]), in order 
% to simplify the code. 
%

nphi=30
ntheta=nphi

[phi,theta,xs,ys,zs]=init_grid_double(nphi,ntheta);

%
% Get the values of the user-defined function 
% on the initial spherical grid (theta in [-pi..pi], see note above)
%
beta = 0;

[xgrid,ygrid,zgrid]=xyz_grid(beta,nphi,xs,ys,ntheta,zs,theta);
fgrid = funuser(xgrid,ygrid,zgrid);

%figure(1)
%scatter3(xgrid,ygrid,zgrid)
%xlabel('x');
%ylabel('y');
%zlabel('z');


%
% Construct the rotated grid. Test an arbitrary ntheta grid.
%
% The algorithm does not assume that the rotated points are
% structured in any way, the user can specify arbitrary number
% and locations of points. In this test, we check a regular
% spherical grid.
%

nphi1=nphi
ntheta1=ntheta/2

[phi1,theta1,xs1,ys1,zs1]=init_grid_single(nphi1,ntheta1);

beta = pi/3;

[xrot,yrot,zrot]=xyz_grid(beta,nphi1,xs1,ys1,ntheta1,zs1,theta1);
frot = funuser(xrot,yrot,zrot);

%figure(2)
%scatter3(xrot,yrot,zrot)
%xlabel('x');
%ylabel('y');
%zlabel('z');

nrot=size(xrot,1)*size(xrot,2)


%
% Get the Fourier coefficients, then, get all rotated grids.
%
tic
fmodes=rotgrid_ref_init(nphi,phi,ntheta,theta,fgrid);
grids=rotgrid_ref_cmpl(nphi,phi,ntheta,theta,fmodes,nrot,xrot,yrot,zrot);
toc


% 
% Finally, check the errors for all rotations
%
errors=zeros(1,nphi);
for j=1:nphi

% Construct the j-th rotated grid directly
  alpha=2*pi*(j-1)/nphi;
  [xrota,yrota,zrota]=xyz_grid_a(beta,nphi1,xs1,ys1,ntheta1,zs1,theta1,alpha);
  frota = funuser(xrota,yrota,zrota);

% Retrieve the rotated grid
  freca=reshape(grids(:,j),nphi1,ntheta1);
 
  errors(j)=norm(frota-freca,2);
end

errors
