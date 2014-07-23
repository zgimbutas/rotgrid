function [phi,theta,xs,ys,zs,ws]=init_grid_lege_double(nphi,ntheta)
%
% Precompute parameters for constructing a double spherical Legendre grid
%
% NOTE: we are using great circles theta = [-pi..pi] in this function
%

phi=2*pi*(0:nphi-1)/nphi;
xs=cos(phi);
ys=sin(phi);
ws=2*pi*ones(1,nphi)/nphi;

% uniformly spaced nodes
%theta=((1:ntheta)-0.5)*2*pi/ntheta;
%zs=cos(theta);

% legendre nodes
ts=grule(ntheta/2); ts=fliplr(ts);

zs=[ts fliplr(ts)];
theta=[acos(ts) 2*pi-fliplr(acos(ts))];

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
