function [x,w,u,v]=chebexps(n)
%CHEBEXPS Chebychev quadrature on the interval [-1,1].
%
%  This subroutine constructs the chebychev nodes on the interval [-1,1],
%  and the weights for the corresponding order n quadrature. it also
%  constructs the matrix v converting the coefficients of a chebychev
%  expansion into its values at the n chebychev nodes, and its inverse u,
%  converting the values of a function at n chebychev nodes into the
%  coefficients of the corresponding chebychev series.  No attempt has
%  been made to make this code efficient, but its speed is normally
%  sufficient.
%
%  X = CHEBEXPS(N) returns the nodes of the Chebychev quadrature rule
%
%  [X,W] = CHEBEXPS(N) returns the nodes and weigths of the Chebychev 
%  quadrature rule
%
%  [X,W,U,V] = CHEBEXPS(N) returns the nodes, weigths and conversion matrices
%  of the Chebychev quadrature rule
%
%  Input parameters:
%
%  N - the number of chebychev nodes and weights to be generated.
%  
%  Output parameters:
%
%  X - the order n chebychev nodes.
%  W - the corresponding quadrature weights.
%  U - the n*n matrix converting the  values at of a polynomial of order
%         n-1 at n chebychev nodes into the coefficients of its 
%         chebychev expansion.
%  V - the n*n matrix converting the coefficients
%         of an n-term chebychev expansion into its values at
%         n chebychev nodes (note that v is the inverse of u).
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

mex_id_ = 'chebexps(i int[x], i int[x], io double[x], io double[xx], io double[xx], io double[x])';
[x, u, v, w] = rotgrid_r2014a(mex_id_, itype, n, x, u, v, w, 1, 1, n, n, n, n, n, n);


