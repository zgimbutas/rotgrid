function [amatr,ts]=legematrin(n,m,xs)
%LEGEMATRIN Interpolation matrix for Gauss-Legendre nodes.
%
%  This subroutine constructs the matrix interpolating
%  functions from the n-point Gaussian grid on the interval [-1,1]
%  to an arbitrary m-point grid (the nodes of the latter are user-provided).
%
%  AMATR = LEGEMATRIN(N,M,XS) 
%
%  Input parameters:
%
%  N - the number of gaussian nodes and weights to be generated.
%  XS - the points at which the function is to be interpolated.
%  
%  Output parameters:
%
%  TS - the order n gaussian nodes.
%

ts=zeros(n,1);
amatr=zeros(m,n);

w=zeros(2*n*n+n+100,1);

mex_id_ = 'lematrin(i int[x], i int[x], i double[], io double[], io double[], io double[])';
[amatr, ts, w] = rotgrid_r2014a(mex_id_, n, m, xs, amatr, ts, w, 1, 1);



