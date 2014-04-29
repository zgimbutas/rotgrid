function [xgrid,ygrid,zgrid]=xyz_grid(beta,nphi,xs,ys,ntheta,zs,theta)
%XYZ_GRID
%
% Construct a spherical grid and rotate it by angle beta around y-axis
%
% Input parameters:
%
% beta - the angle of rotation around y-axis
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

xgrid = zeros(nphi,ntheta);
ygrid = zeros(nphi,ntheta);
zgrid = zeros(nphi,ntheta);

costheta=cos(theta);
sintheta=sin(theta);

if( 1 == 2 ),
for i=1:nphi,
  for j=1:ntheta,
%    z=zs(j);
%    x=xs(i)*sqrt(1-z**2);
%    y=ys(i)*sqrt(1-z**2);
    z=costheta(j);
    x=+xs(i)*sintheta(j);
    y=+ys(i)*sintheta(j);
%  counter-clockwise rotation: positive angle
%    xgrid(i,j) = x*cosbeta - z*sinbeta;
%    ygrid(i,j) = y;
%    zgrid(i,j) = x*sinbeta + z*cosbeta;
%  clockwise rotation: negative angle
    xgrid(i,j) = x*cosbeta + z*sinbeta;
    ygrid(i,j) = y;
    zgrid(i,j) =-x*sinbeta + z*cosbeta;
  end
end
end

if( 2 == 2 ),
x=kron(sintheta,xs');
y=kron(sintheta,ys');
z=kron(costheta,ones(nphi,1));

xgrid = x*cosbeta - z*sinbeta;
ygrid = y;
zgrid = x*sinbeta + z*cosbeta;

%  counter-clockwise rotation: positive angle
%    xgrid = x*cosbeta - z*sinbeta;
%    ygrid = y;
%    zgrid = x*sinbeta + z*cosbeta;
%  clockwise rotation: negative angle
    xgrid = x*cosbeta + z*sinbeta;
    ygrid = y;
    zgrid =-x*sinbeta + z*cosbeta;

end

