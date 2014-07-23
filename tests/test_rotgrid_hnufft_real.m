%  Testing routines for fast gridding algorithm
%
%  Fast rotation of the user-defined grid (xrot,yrot,zrot).
%  This function constructs a set of grids, obtained by rotating 
%  the original grid by nphi uniformly spaced angles phi.


addpath '../nufft1d';
addpath '../nufft2d';


p = 12*4;
ntheta0=p+1;

%
% Construct the initial spherical grid, 
% it will be used to get Fourier coefficients for a sphere.
%
% NOTE: we are using big circles theta = [-pi..pi] 
% (and not the meridians theta = [0..pi]), in order 
% to simplify the code. 
%

ntheta=ntheta0*2;
nphi=ntheta0*2;

nphi=fftnext235(nphi);

[phi,theta,xs,ys,zs]=init_grid_cheb_double(nphi,ntheta);

%
% Get the values of the user-defined function 
% on the initial spherical grid (theta in [-pi..pi], see note above)
%
beta = 0;

[xgrid,ygrid,zgrid]=xyz_grid(beta,nphi,xs,ys,ntheta,zs,theta);
fgrid = funuser_real(xgrid,ygrid,zgrid);

%figure(1)
%scatter3(xgrid,ygrid,zgrid)
%xlabel('x');
%ylabel('y');
%zlabel('z');


%
% Construct the rotated grid. Test an arbitrary ntheta grid.
%

disp('')
disp('TEST #1, rotgrid_real')

% The algorithm does not assume that the rotated points are
% structured in any way, the user can specify arbitrary number
% and locations of points. 


nphi1=nphi;
ntheta1=ntheta/2;

ngrid = nphi1*ntheta1

% randomly distributed points on a unit sphere
theta=rand(ngrid,1)*pi;
phi=rand(ngrid,1)*2*pi;
r=1;
xgrid=r*cos(phi).*sin(theta);
ygrid=r*sin(phi).*sin(theta);
zgrid=r*cos(theta);

%
% Get the Fourier coefficients, then, get all rotated grids.
%
tic
fmodes=rotgrid_init_dcheb(nphi,phi,ntheta,theta,fgrid);
toc
tic
grids=rotgrid_real(nphi,phi,ntheta,theta,fmodes,ngrid,xgrid,ygrid,zgrid);
toc


% 
% Finally, check the errors for all rotations
%
errors=zeros(1,nphi);
for j=1:nphi

% Construct the j-th rotated grid directly
  alpha=2*pi*(j-1)/nphi;
  [xrota,yrota,zrota]=rotgrid_direct_ba(beta,alpha,xgrid,ygrid,zgrid);
  frota = funuser_real(xrota,yrota,zrota);

% Retrieve the rotated grid
  freca=grids(:,j);
 
  errors(j)=norm(frota-freca,2);
end

max(abs(errors))


disp('')
disp('TEST #2, rotgrid_real_opt')

%
% In this test, we check a regular spherical Legendre grid.
%

nphi1=nphi
ntheta1=ntheta/2

[phi1,theta1,xs1,ys1,zs1]=init_grid_lege_single(nphi1,ntheta1);

beta = pi/3;

[xgrid,ygrid,zgrid]=xyz_grid(beta,nphi1,xs1,ys1,ntheta1,zs1,theta1);

%figure(2)
%scatter3(xgrid,ygrid,zgrid)
%xlabel('x');
%ylabel('y');
%zlabel('z');

ngrid=size(xgrid,1)*size(xgrid,2)

%
% Get the Fourier coefficients, then, get all rotated grids.
%
tic
fmodes=rotgrid_init_dcheb(nphi,phi,ntheta,theta,fgrid);
toc
tic
grids=rotgrid_real_opt(nphi,phi,ntheta,theta,fmodes,ngrid,xgrid,ygrid,zgrid);
toc


% 
% Finally, check the errors for all rotations
%
errors=zeros(1,nphi);
for j=1:nphi

% Construct the j-th rotated grid directly
  alpha=2*pi*(j-1)/nphi;
  [xrota,yrota,zrota]=xyz_grid_a(beta,nphi1,xs1,ys1,ntheta1,zs1,theta1,alpha);
  frota = funuser_real(xrota,yrota,zrota);

% Retrieve the rotated grid
  freca=reshape(grids(:,j),nphi1,ntheta1);
 
  errors(j)=norm(frota-freca,2);
end

max(abs(errors))


disp('')
disp('TEST #3, rotgrid_real_opt_vec')

%
% In this test, we check a regular spherical Legendre grid.
% Merge grids for multiple latitudes.
%
ngrid= nphi1*ntheta1
npoints = ngrid*ntheta1

xgrid= zeros(nphi1,ntheta1,ntheta1);
ygrid= zeros(nphi1,ntheta1,ntheta1);
zgrid= zeros(nphi1,ntheta1,ntheta1);
for i=1:ntheta1
[xrot1,yrot1,zrot1]=xyz_grid(theta1(i),nphi1,xs1,ys1,ntheta1,zs1,theta1);
xgrid(:,:,i)=xrot1(:,:);
ygrid(:,:,i)=yrot1(:,:);
zgrid(:,:,i)=zrot1(:,:);
end


%
% Get the Fourier coefficients, then, get all rotated grids.
%
tic
fmodes=rotgrid_init_dcheb(nphi,phi,ntheta,theta,fgrid);
toc
tic
grids1=rotgrid_real_opt_vec(nphi,phi,ntheta,theta,fmodes,ngrid,ntheta1,xgrid,ygrid,zgrid);
toc
grids1=reshape(grids1,ngrid,nphi,ntheta1);


% 
% Finally, check the errors for all rotations
%
errors=zeros(1,nphi);
k=3;
for j=1:nphi

% Construct the j-th rotated grid directly
  alpha=2*pi*(j-1)/nphi;
  [xrota,yrota,zrota]=...
      xyz_grid_a(theta1(k),nphi1,xs1,ys1,ntheta1,zs1,theta1,alpha);
  frota = funuser_real(xrota,yrota,zrota);

% Retrieve the rotated grid
  freca=reshape(grids1(:,j,k),nphi1,ntheta1);
 
  errors(j)=norm(frota-freca,2);
end

max(abs(errors))
