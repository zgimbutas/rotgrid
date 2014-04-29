function [x,w]=sphere_singular(nmax,itype)
%SHPERE_SINGULAR: Quadratures for smooth and 1/r singularities on the sphere

if( itype == 1 ),
%  weights for smooth functions on the sphere
% 
  n = nmax+1;
  [x,w]=legeexps(n);
  x=-x;
end

if( itype == 2 ),
%  weights for smooth functions on the sphere
% 
  n = 2*nmax-1;
  [x,w]=chebexps(n);
  x=-x;
end

if( itype == 3 ),
%  weights for convolution with 1/r on the sphere
% 
  n = nmax+1;
  [x,w]=legeexps(n);
  x=-x;
  for i=1:n
    pols=legepols(x(i),n-1);
    d=sum(pols);
    w(i)=w(i)*d;
  end
end

if( itype == 4 ),
%  weights for convolution with 1/r on the sphere
% 
  n = 2*nmax-1;
  [x,w]=chebexps(n);
  x=-x;
  for i=1:n
    pols=legepols(x(i),n-1);
    d=sum(pols);
    w(i)=w(i)*d;
  end
end


if( itype == 5 ),
%  construct the Jacobi nodes and weights on the interval [-1,1]
%
 n = nmax+1;
 ikind=5;
 kpts=0; endpts=[-1;1];
 alpha=-0.5;  beta=0.0;
 [x,w] = gaussq(ikind,n,alpha,beta,kpts,endpts);
 x=flipud(x);  w=flipud(w);
 w=w*sqrt(2.0)/2;
end


if( itype == 6 ),
%  weights for convolution with 1/r + smooth functions on the sphere
% 
  n = nmax+1;

% construct the generalized Gaussian quadrature
% for P_n(x) + 1/sqrt(x) P_n(x) on interval [0,1]
%
  [x0,w0]=legeexps(n);
  x=((x0+1)/2).^2;
  w=w0.*(x0+1)/2;

% construct the generalized Gaussian quadrature
% for P_n(x) + 1/sqrt(1-x) P_n(x) on interval [-1,1]
% (north pole)
%
  x=-(x*2-1);
  w=w*2;
end

