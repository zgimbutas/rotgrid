function grids=rotgrid_dsr_real(nterms,mpole,nphi,ntheta,nbeta,beta,nrot,rotmat,ctheta,ynms,wsave)
%ROTGRID_FSR_REAL Rotate the real spherical harmonic grids.
%
%  Fast, direct algorithm for rotating real spherical harmonic grids
%  (dimensioned NPHI-by-NTHETA) into new pole locations (beta_j, alpha_k), 
%  where alpha_k = 2*pi * k/nrot, k=0..nrot-1,  beta_j, j=1..nbeta. 
%
%  GRIDS = rotgrid_fsr_real(NTERMS,MPOLE,NPHI,NTHETA,NBETA,BETA,NROT,...
%       ROTMAT,CTHETA,YNMS,WSAVE) rotates 
%  the real spherical harmonics expansion of degree NTERMS 
%  about the z-axis by degree ALPHA_K and about the y-axis by degree BETA.
%  into a collection of new pole locations (beta_j, alpha_k) in spherical 
%  coordinates (theta, phi), where alpha_k = 2*pi * k/nrot, k=0..nrot-1,
%  beta_j, j=1..nbeta.
%
%  The rotated poles form a uniformly spaced grid on lattitude \theta.
%
%  grids - function values on the rotated grids, 
%             NPHI-by-NTHETA-by-NROT-by-NBETA real*8 matrix
%
%      Input parameters:
%
%  nterms - the number of terms in spherical harmonics expansion
%  mpole - the coefficients of spherical harmonics expansion,
%                    complex*16 (0:nterms,0:nterms)
%  nphi - the number of points in latitude discretization (for spherical grid)
%  ntheta - the number of points in meridian discretization (for spherical grid)
%  ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
%  ynms - Legendre functions, real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
%  wsave - fftpack wsave array, complex *16 wsave(4*nphi+15)
%
%  nbeta - the number of points in pole location meridian discretization
%  beta - angles for new pole locations beta_j, j=1..nbeta
%  nrot - angles for new pole locations alpha_k = 2*pi * k/nrot, k=0..nrot-1.
%       
%  rotmat - The rotation operators for directions (beta_j,0). 
% 
%  [rotmat] must be initialized via 
%              a preceding call to rotgrid_fsr_real_init
%  [ctheta, ynms, wsave] must be initialized via 
%      a preceding call to sphtrans_real_lege_init or sphtrans_real_cheb_init
%
%      Output parameters:
%
%  grids - function values on the rotated grids, 
%             NPHI-by-NTHETA-by-NROT-by-NBETA real*8 matrix
%
%  Our definition of complex spherical harmonics is
%
%  Ynm(theta,phi)= sqrt( 2n+1) sqrt((n-m)!/(n+m)!) 
%                  Pnm(cos theta) e^(im phi), 
%  Yn,-m(theta,phi) = sqrt( 2n+1) sqrt((n-m)!/(n+m)!) 
%                  Pnm(cos theta) e^(-im phi),   for m >= 0.
%       
%  Note that we do not include the Condon-Shortley phase (-1)^m, if m<0.
%

grids=zeros(nphi*ntheta,nrot*nbeta);

mex_id_ = 'rotgrid_dsr_real(i int[], i dcomplex[], i int[], i int[], i int[], i double[], i int[], io double[], i double[], i double[], i double[], i dcomplex[])';
[grids] = rotgrid_r2014a(mex_id_, nterms, mpole, nphi, ntheta, nbeta, beta, nrot, grids, rotmat, ctheta, ynms, wsave);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


