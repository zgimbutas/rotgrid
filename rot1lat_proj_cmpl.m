function mpout=rot1lat_proj_cmpl(nterms,mpole,beta,nrot)
%ROT1LAT_PROJ_CMPL Rotate the complex spherical harmonics expansion.
%
%  MPOUT = rot1lat_proj_cmpl(NTERMS,MPOLE,BETA,NROT) rotates 
%  the complex spherical harmonics expansion of degree NTERMS 
%  about the z-axis by degree ALPHA_K and about the y-axis by degree BETA.
%  into a collection of new pole locations (beta, alpha_k) in spherical 
%  coordinates (theta, phi), where alpha_k = 2*pi * k/nrot, k=0..nrot-1.
%
%  The rotated poles form a uniformly spaced grid on lattitude \theta.
%
%  MPOLE is (NTERMS+1)-by-(2*NTERMS+1) complex matrix.
%  MPOUT is (NTERMS+1)-by-(2*NTERMS+1)-by-NROT complex matrix.
%
%
%  Fast and stable algorithm for rotating complex spherical harmonic expansions
%  into locations (beta, alpha_k), where alpha_k = 2*pi * k/nrot, k=0..nrot-1.
%
%  The method is based on computing the induced potential and
%  its theta-derivative on the rotated equator
%  for each order (first index). The coefficients of  the rotated
%  expansion can then be obtained by FFT and projection.
%
%  There is some loss in speed over using recurrence relations 
%  but it is stable to all orders whereas the recurrence schemes 
%  are not.
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

if( nargin == 0 ), mpout=rot1lat_proj_cmpl_test(); return; end;

m1=nterms;
m2=nterms;

%%%mpout = zeros(nterms+1,2*nterms+1,nrot)+1i*zeros(nterms+1,2*nterms+1,nrot);
mpout = zeros((nterms+1)*(2*nterms+1),nrot)*(1+1i);

lmp = nterms;
lmpn = nterms;

mex_id_ = 'rotviaproj3f90(i double[x], i int[x], i int[x], i int[x], i int[x], i dcomplex[], i int[x], io dcomplex[], i int[x])';
[mpout] = rotgrid_r2014a(mex_id_, beta, nrot, nterms, m1, m2, mpole, lmp, mpout, lmpn, 1, 1, 1, 1, 1, 1, 1);

mpout=reshape(mpout,nterms+1,2*nterms+1,nrot);

function mpout=rot1lat_proj_cmpl_test()

nterms = 2;
beta = pi/2;
nrot = 2;

mpole = zeros(nterms+1,2*nterms+1)+1i*zeros(nterms+1,2*nterms+1);

j=nterms+1;
mpole(1,j) = 1;
mpole(2,j+(-1:1)) = 1; 
mpole(3,j+(-2:2)) = 1; 

mpole;

mpout=rot1lat_proj_cmpl(nterms,mpole,beta,nrot);
mpout=reshape(mpout,nterms+1,2*nterms+1,nrot);


