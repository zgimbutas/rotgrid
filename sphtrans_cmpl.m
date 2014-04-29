function fgrid=sphtrans_cmpl(nterms,mpole,nphi,ntheta,ctheta,ynms,wsave)
%
%  Complex valued O(p^3) spherical transform on a spherical grid. 
%  Assumes a symmetric grid, only half of Legendre functions are stored.
%
%  Input parameters:
%
%  nterms - the number of terms in spherical harmonics expansion
%  mpole - the coefficients of spherical harmonics expansion,
%               complex*16 (0:nterms,-nterms:nterms)
%  nphi - the number of points in latitude discretization
%  ntheta - the number of points in meridian discretization
%  ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
%  ynms - Legendre functions, real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
%  wsave - fftpack wsave array, complex *16 wsave(4*nphi+15)
%       
%  Output parameters:
%
%  fgrid - function values on the grid, NPHI-by-NTHETA complex*16 matrix
%
fgrid = zeros(nphi,ntheta)+1i*zeros(nphi,ntheta);

mex_id_ = 'sphtrans_cmpl(i int[], i dcomplex[], i int[], i int[], io dcomplex[], i double[], i double[], i dcomplex[])';
[fgrid] = rotgrid_r2014a(mex_id_, nterms, mpole, nphi, ntheta, fgrid, ctheta, ynms, wsave);


