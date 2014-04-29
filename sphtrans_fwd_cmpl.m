function mpole=sphtrans_fwd_cmpl(nterms,nphi,ntheta,fgrid,ctheta,whts,ynms,wsave)
%
%  Complex valued O(p^3) forward spherical transform on a spherical grid. 
%  Assumes a symmetric grid, only half of Legendre functions are stored.
%
%  Input parameters:
%
%  nterms - the number of terms in spherical harmonics expansion
%  nphi - the number of points in latitude discretization
%  ntheta - the number of points in meridian discretization
%  fgrid - function values on the grid, NPHI-by-NTHETA complex*16 matrix
%  ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
%  whts - weights of meridian discretization angles, real*8 ntheta
%  ynms - Legendre functions, real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
%  wsave - fftpack wsave array, complex *16 wsave(4*nphi+15)
%       
%  Output parameters:
%
%  mpole - the coefficients of spherical harmonics expansion,
%               complex*16 (0:nterms,-nterms:nterms)
%
mpole = zeros(nterms+1,2*nterms+1)+1i*zeros(nterms+1,2*nterms+1);

mex_id_ = 'sphtrans_fwd_cmpl(i int[], io dcomplex[], i int[], i int[], i dcomplex[], i double[], i double[], i double[], i dcomplex[])';
[mpole] = rotgrid_r2014a(mex_id_, nterms, mpole, nphi, ntheta, fgrid, ctheta, whts, ynms, wsave);


