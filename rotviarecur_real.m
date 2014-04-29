function mpout=rotviarecur_real(nterms,mpole,beta,m1,m2)
%ROTVIARECUR_REAL Rotate the real spherical harmonics about the y-axis. 
%
%  MPOUT = rotviarecur_real(NTERMS,MPOLE,BETA) rotates the real spherical 
%  harmonics expansion of degree NTERMS about the y-axis by degree BETA. 
%
%  MPOUT = rotviarecur_real(NTERMS,MPOLE,BETA,M1,M2) rotates the real spherical 
%  harmonics expansion of degree NTERMS about the y-axis by degree BETA. 
%  Only modes up to degree M1 in the input MPOLE expansion and modes up 
%  to degree M2 in the output expansion MPOUT will be used. 
%  
%  Both MPOUT and MPOLE are (NTERMS+1)-by-(NTERMS+1) complex matrices.
%
%
%  Fast, recursive algorithm for applying rotation operator about
%  the y-axis determined by angle beta. It is good for NTERMS up to 100 or so.
%
%  Our definition of real spherical harmonics is
%
%  Ynm(theta,phi)= sqrt( 2n+1) sqrt((n-m)!/(n+m)!) 
%                  Pnm(cos theta) e^(im phi), 
%  Yn,-m(theta,phi) = sqrt( 2n+1) sqrt((n-m)!/(n+m)!) 
%                  Pnm(cos theta) e^(-im phi),   for m >= 0.
%       
%  Note that we do not include the Condon-Shortley phase (-1)^m, if m<0.
%

if( nargin == 0 ), mpout=rotviarecur_real_test(); return; end;

if( nargin == 3 ),
  m1=nterms;
  m2=nterms;
end;

mpout = zeros(nterms+1,nterms+1)+1i*zeros(nterms+1,nterms+1);

lmp = nterms;
lmpn = nterms;

mex_id_ = 'rotviarecur3f90_real(i double[x], i int[x], i int[x], i int[x], io dcomplex[], i int[x], io dcomplex[], i int[x])';
[mpole, mpout] = rotgrid_r2014a(mex_id_, beta, nterms, m1, m2, mpole, lmp, mpout, lmpn, 1, 1, 1, 1, 1, 1);


function mpout=rotviarecur_real_test()

nterms = 2;
beta = pi/2;
nrot = 2;

mpole = zeros(nterms+1,nterms+1)+1i*zeros(nterms+1,nterms+1);

j=1;
mpole(1,j) = 1;
mpole(2,j+(0:1)) = 1; 
mpole(3,j+(0:2)) = 1; 

mpout = rotviarecur_real(nterms,mpole,beta);




