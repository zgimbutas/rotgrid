function [ctheta,whts,ynms,wsave]=sphtrans_cmpl_lege_init(nterms,nphi,ntheta)
%
%  Precompute parameters and tables for O(p^3) spherical transforms
%  on Legendre spherical grid. Assumes a symmetric grid, only half of
%  Legendre functions are stored.
%
%  Input parameters:
%
%  nterms - the number of terms in spherical harmonics expansion
%  nphi - the number of points in latitude discretization
%  ntheta - the number of points in meridian discretization
%       
%  Output parameters:
%       
%  ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
%  whts - weights of meridian discretization angles, real*8 ntheta
%  ynms - Legendre functions, real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
%  wsave - initialized fftpack wsave array, complex *16 wsave(4*nphi+15)
%
ctheta = zeros(ntheta,1);
whts = zeros(ntheta,1);
ynms = zeros((nterms+1)^2,ntheta/2+1);
wsave = zeros(4*nphi+15,1)+1i*zeros(4*nphi+15,1);

mex_id_ = 'sphtrans_cmpl_lege_init(i int[], i int[], i int[], io double[], io double[], io double[], io dcomplex[])';
[ctheta, whts, ynms, wsave] = rotgrid_r2014a(mex_id_, nterms, nphi, ntheta, ctheta, whts, ynms, wsave);


