function [phi,theta,xs,ys,zs,ws]=init_grid_double(nphi,ntheta)
%
% Precompute parameters for constructing a double spherical Chebychev grid
%
% NOTE: we are using great circles theta = [-pi..pi] in this function
%

phi=2*pi*(0:nphi-1)/nphi;
xs=cos(phi);
ys=sin(phi);
ws=2*pi*ones(1,nphi)/nphi;

theta=((1:ntheta)-0.5)*2*pi/ntheta;
zs=cos(theta);

%
% Unwrap phases to interval [-pi,pi], if needed
%
if( 1 == 2 ),
for j=1:ntheta,
  if( theta(j) > pi ), theta(j)=theta(j)-2*pi; end
end
end

if( 2 == 2 ),
j=find(theta > pi);
theta(j)=theta(j)-2*pi;
end
