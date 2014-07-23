function fgrid=sphtrans_real(nterms,mpole,nphi,ntheta,ctheta,ynms,wsave)
%
%  Real valued O(p^3) spherical transform on a spherical grid. 
%  Assumes a symmetric grid, only half of Legendre functions are stored.
%  No aliasing (nphi >= 2*nterms+1)
%
%  Input parameters:
%
%  nterms - the number of terms in spherical harmonics expansion
%  mpole - the coefficients of spherical harmonics expansion,
%               complex*16 (0:nterms,0:nterms)
%  nphi - the number of points in latitude discretization
%  ntheta - the number of points in meridian discretization
%  ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
%  ynms - Legendre functions, real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
%  wsave - fftpack wsave array, complex *16 wsave(4*nphi+15)
%       
%  Output parameters:
%
%  fgrid - function values on the grid, NPHI-by-NTHETA real*8 matrix
%
fgrid = zeros(nphi,ntheta);

mex_id_ = 'sphtrans_real(i int[], i dcomplex[], i int[], i int[], io double[], i double[], i double[], i dcomplex[])';
[fgrid] = rotgrid_r2014a(mex_id_, nterms, mpole, nphi, ntheta, fgrid, ctheta, ynms, wsave);



