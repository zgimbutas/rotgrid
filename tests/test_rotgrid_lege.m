%  Testing routines for fast gridding algorithm
%
%  Fast rotation of the user-defined grid (xrot,yrot,zrot).
%  This function constructs a set of grids, obtained by rotating 
%  the original grid by nphi uniformly spaced angles phi.

%  Test Legendre grids.

addpath '../matlab';
addpath '../nufft1d';
addpath '../nufft2d';


%
% Construct the initial double spherical Chebychev grid, 
% it will be used to get Fourier coefficients for a sphere.
%
% NOTE: we are using great circles theta = [-pi..pi] 
% (and not the meridians theta = [0..pi]), in order 
% to simplify the code. 
%

nphi=8
ntheta=nphi
 
%
% Get the values of the user-defined function 
% on the initial spherical grid (theta in [-pi..pi], see note above)
%
beta = 0;

if( 2 == 2 ),
% double grid
[phi,theta,xs,ys,zs]=init_grid_lege_double(nphi,ntheta);

[xgrid,ygrid,zgrid]=xyz_grid(beta,nphi,xs,ys,ntheta,zs,theta);
fgrid = funuser(xgrid,ygrid,zgrid);

end

if( 1 == 2 ),
% single grid
ntheta0=ntheta/2
[phi0,theta0,xs0,ys0,zs0]=init_grid_lege_single(nphi,ntheta0);

[xgrid,ygrid,zgrid]=xyz_grid(beta,nphi,xs0,ys0,ntheta0,zs0,theta0);
fgrid = funuser(xgrid,ygrid,zgrid);

a = fgrid(1:nphi/2,:);
b = fgrid(nphi/2+1:nphi,:);
fgrid = [[a; b] [fliplr(b); fliplr(a)]]

end

% Convert the spherical Legendre grid to the spherical Chebychev grid
[amatr,xc,xg,c,g] = legecheb_double(ntheta);
%%%fgrid = ((amatr * fgrid.').');
fgrid = fgrid * amatr.';
[phi,theta,xs,ys,zs]=init_grid_cheb_double(nphi,ntheta);


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
% spherical Legendre grid.
%

nphi1=nphi
ntheta1=ntheta/2

[phi1,theta1,xs1,ys1,zs1]=init_grid_lege_single(nphi1,ntheta1);

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
fmodes=rotgrid_init_dcheb(nphi,phi,ntheta,theta,fgrid);
grids=rotgrid_cmpl(nphi,phi,ntheta,theta,fmodes,nrot,xrot,yrot,zrot);
toc


% 
% Finnaly, check the errors for all rotations
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
