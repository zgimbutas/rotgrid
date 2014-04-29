function [a,xc,xg,c,g]=legecheb_double(n)
%LEGECHEB_DOUBLE Interpolation matrix for double Legendre grids.
%
% This function constructs the interpolation matrix 
% converting the function values at spherical Legendre double grid 
% into the function values the spherical Chebychev double grid.
%
%  [A,XC,XG,C,G]=LEGECHEB_DOUBLE(N)
%
%  Input parameters:
%
%  N - the number of nodes in double Legendre grid to be generated, 
%      must be even for double grids.
%  
%  Output parameters:
%
%  A - the n*n interpolation matrix converting the function values 
%    at double Legendre grid into the function values at double Chebychev grid.
%  XC - uniformly spaced Chebychev nodes for double grid.
%  XG - Legendre nodes for double grid.
%  C - the n*n matrix converting the coefficients of Fourier expansion
%      into the function values at double Chebychev grid.
%  G - the n*n matrix converting the coefficients of Fourier expansion
%      into the function values at double Legendre grid.
%

xc = 2*pi*((1:n)-0.5)/n;
c = exp(1i*kron(-n/2:(n-1)/2,xc'));

[ts]=grule(n/2); ts=fliplr(ts);
xg=[acos(ts) 2*pi-fliplr(acos(ts))];
g = exp(1i*kron(-n/2:(n-1)/2,xg'));

a = c*inv(g);
