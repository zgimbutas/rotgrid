function [xrot,yrot,zrot]=rotgrid_direct_ba(beta,alpha,xgrid,ygrid,zgrid)
%ROTGRID_DIRECT_BA
%
% Rotate a spherical grid by angle beta around y-axis,
% then, rotate it by angle alpha around z-axis
%
% Input parameters:
%
% beta - the angle of rotation around y-axis
% alpha - the angle of rotation around z-axis
% xgrid - the x-coordinates of the constructed grid
% ygrid - the y-coordinates of the constructed grid
% zgrid - the z-coordinates of the constructed grid
%
% Output parameters:
%
% xrot - the x-coordinates of the rotated grid
% yrot - the y-coordinates of the rotated grid
% zrot - the z-coordinates of the rotated grid
%

cosbeta=cos(beta);
sinbeta=sin(beta);

sinalpha=sin(alpha);
cosalpha=cos(alpha);

%  counter-clockwise rotation: positive angle
%    x1 = xgrid*cosbeta - zgrid*sinbeta;
%    y1 = ygrid;
%    z1 = xgrid*sinbeta + zgrid*cosbeta;

%  clockwise rotation: negative angle
    x1 = xgrid*cosbeta + zgrid*sinbeta;
    y1 = ygrid;
    z1 =-xgrid*sinbeta + zgrid*cosbeta;

xrot = x1*cosalpha - y1*sinalpha;
yrot = x1*sinalpha + y1*cosalpha;
zrot = z1;


