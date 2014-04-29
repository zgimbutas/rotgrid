function mpout=rot1lat_fsr_real(nterms,mpole,beta,nrot,rotmat,ldc)
%ROT1LAT_FSR_REAL Rotate the real spherical harmonics expansion.
%
%  MPOUT = rot1lat_fsr_real(NTERMS,MPOLE,BETA,NROT,ROTMAT,LDC) rotates 
%  the real spherical harmonics expansion of degree NTERMS 
%  about the z-axis by degree ALPHA_K and about the y-axis by degree BETA.
%  into a collection of new pole locations (beta, alpha_k) in spherical 
%  coordinates (theta, phi), where alpha_k = 2*pi * k/nrot, k=0..nrot-1.
%
%  The rotated poles form a uniformly spaced grid on lattitude \theta.
%
%  MPOLE is (NTERMS+1)-by-(NTERMS+1) complex matrix.
%  MPOUT is (NTERMS+1)-by-(NTERMS+1)-by-NROT complex matrix.
%
%
%  Fast, FFT-based algorithm for rotating real spherical harmonic expansions
%  into locations (beta, alpha_k), where alpha_k = 2*pi * k/nrot, k=0..nrot-1.
%
%  rotmat - The rotation operator for direction (beta,0). It must be initialized
%     via a preceding call to real(rotmat_proj_init(ldc,beta)), 
%     Note, that this function expects real-valued rotation matrices, please
%     take the real part of rotmat after the call to rotmat_proj_init, see
%     the testing code below.
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

if( nargin == 0 ), mpout=rot1lat_fsr_real_test(); return; end;

m1=nterms;
m2=nterms;

%%%mpout = zeros(nterms+1,nterms+1,nrot)+1i*zeros(nterms+1,nterms+1,nrot);
mpout = zeros((nterms+1)*(nterms+1),nrot)*(1+1i);

lmp = nterms;
lmpn = nterms;

mex_id_ = 'rot1lat_wfft_real(i double[x], i int[x], i int[x], i int[x], i int[x], io dcomplex[], i int[x], io dcomplex[], i int[x], i double[], i int[x])';
[mpole, mpout] = rotgrid_r2014a(mex_id_, beta, nrot, nterms, m1, m2, mpole, lmp, mpout, lmpn, rotmat, ldc, 1, 1, 1, 1, 1, 1, 1, 1);

mpout=reshape(mpout,nterms+1,nterms+1,nrot);


function mpout=rot1lat_fsr_real_test()

nterms = 2;
beta = pi/2;
nrot = 2;

mpole = zeros(nterms+1,nterms+1)+1i*zeros(nterms+1,nterms+1);

j=1;
mpole(1,j) = 1;
mpole(2,j+(0:1)) = 1; 
mpole(3,j+(0:2)) = 1; 

mpole;

ldc = nterms;
rotmat = real(rotmat_proj_init(ldc,beta));

mpout=rot1lat_fsr_real(nterms,mpole,beta,nrot,rotmat,ldc);
mpout=reshape(mpout,nterms+1,nterms+1,nrot);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


