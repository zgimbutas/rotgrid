function [xgrid,ygrid,zgrid]=xyz_grid_a(beta,nphi,xs,ys,ntheta,zs,theta,alpha)
%XYZ_GRID_A
%
% Construct a spherical grid and rotate it by angle beta around y-axis,
% then, rotate it by angle alpha around z-axis
%
% Input parameters:
%
% beta - the angle of rotation around y-axis
% alpha - the angle of rotation around z-axis
% nphi, ntheta - the number of discretization points in phi and theta
% phi,theta,xs,ys,zs - must be constructed via a preceding call to 
%                      either init_grid_single or init_grid_double
%
% Output parameters:
%
% xgrid - the x-coordinates of the constructed grid
% ygrid - the y-coordinates of the constructed grid
% zgrid - the z-coordinates of the constructed grid
%

cosbeta=cos(beta);
sinbeta=sin(beta);

sinalpha=sin(alpha);
cosalpha=cos(alpha);

xgrid = zeros(nphi,ntheta);
ygrid = zeros(nphi,ntheta);
zgrid = zeros(nphi,ntheta);

costheta=cos(theta);
sintheta=sin(theta);

if( 1 == 2 ),
for i=1:nphi,
  for j=1:ntheta,
    z=costheta(j);
    x=+xs(i)*sintheta(j);
    y=+ys(i)*sintheta(j);
%  counter-clockwise rotation: positive angle
%    x1 = x*cosbeta - z*sinbeta;
%    y1 = y;
%    z1 = x*sinbeta + z*cosbeta;
%  clockwise rotation: negative angle
    x1 = x*cosbeta + z*sinbeta;
    y1 = y;
    z1 =-x*sinbeta + z*cosbeta;
    xgrid(i,j) = x1*cosalpha - y1*sinalpha;
    ygrid(i,j) = x1*sinalpha + y1*cosalpha;
    zgrid(i,j) = z1;
  end
end
end

if( 2 == 2 ),
x=kron(sintheta,xs');
y=kron(sintheta,ys');
z=kron(costheta,ones(nphi,1));

%  counter-clockwise rotation: positive angle
%    x1 = x*cosbeta - z*sinbeta;
%    y1 = y;
%    z1 = x*sinbeta + z*cosbeta;
%  clockwise rotation: negative angle
    x1 = x*cosbeta + z*sinbeta;
    y1 = y;
    z1 =-x*sinbeta + z*cosbeta;

xgrid = x1*cosalpha - y1*sinalpha;
ygrid = x1*sinalpha + y1*cosalpha;
zgrid = z1;
end

