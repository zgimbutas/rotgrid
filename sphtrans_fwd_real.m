function mpole=sphtrans_fwd_real(nterms,nphi,ntheta,fgrid,ctheta,whts,ynms,wsave)
%
%  Real valued O(p^3) forward spherical transform on a spherical grid. 
%
%  Input parameters:
%
%  nterms - the number of terms in spherical harmonics expansion
%  nphi - the number of points in latitude discretization
%  ntheta - the number of points in meridian discretization
%  fgrid - function values on the grid, NPHI-by-NTHETA real*8 matrix
%  ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
%  whts - weights of meridian discretization angles, real*8 ntheta
%  ynms - Legendre functions, real *8 ynms(0:nterms,0:nterms,ntheta)
%  wsave - fftpack wsave array, complex *16 wsave(4*nphi+15)
%       
%  Output parameters:
%
%  mpole - the coefficients of spherical harmonics expansion,
%               complex*16 (0:nterms,0:nterms)
%
mpole = zeros(nterms+1,nterms+1)+1i*zeros(nterms+1,nterms+1);

mex_id_ = 'sphtrans_fwd_real(i int[], io dcomplex[], i int[], i int[], i double[], i double[], i double[], i double[], i dcomplex[])';
[mpole] = rotgrid_r2014a(mex_id_, nterms, mpole, nphi, ntheta, fgrid, ctheta, whts, ynms, wsave);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


