function [x,w,u,v]=legeexps(n)
%LEGEEXPS Gauss-Legendre quadrature on the interval [-1,1].
%
%  This subroutine constructs the gaussian nodes on the interval [-1,1],
%  and the weights for the corresponding order n quadrature. it also
%  constructs the matrix v converting the coefficients of a legendre
%  expansion into its values at the n gaussian nodes, and its inverse u,
%  converting the values of a function at n gaussian nodes into the
%  coefficients of the corresponding legendre series.  No attempt has
%  been made to make this code efficient, but its speed is normally
%  sufficient.
%
%  X = LEGEEXPS(N) returns the nodes of the Gauss-Legendre quadrature rule
%
%  [X,W] = LEGEEXPS(N) returns the nodes and weigths of the Gauss-Legendre 
%  quadrature rule
%
%  [X,W,U,V] = LEGEEXPS(N) returns the nodes, weigths and conversion matrices
%  of the Gauss-Legendre quadrature rule
%
%  Input parameters:
%
%  N - the number of gaussian nodes and weights to be generated.
%  
%  Output parameters:
%
%  X - the order n gaussian nodes.
%  W - the corresponding quadrature weights.
%  U - the n*n matrix converting the  values at of a polynomial of order
%         n-1 at n legendre nodes into the coefficients of its 
%         legendre expansion.
%  V - the n*n matrix converting the coefficients
%         of an n-term legendre expansion into its values at
%         n legendre nodes (note that v is the inverse of u).
%

x=zeros(1,n);
w=zeros(1,n);
u=zeros(n,n);
v=zeros(n,n);

itype=0;
if( nargout == 1 ) itype=0; end;
if( nargout == 2 ) itype=1; end;
if( nargout == 3 ) itype=2; end;
if( nargout == 4 ) itype=2; end;

mex_id_ = 'legeexps(i int[x], i int[x], i double[x], i double[xx], i double[xx], i double[x])';
rotgrid_r2014a(mex_id_, itype, n, x, u, v, w, 1, 1, n, n, n, n, n, n);


