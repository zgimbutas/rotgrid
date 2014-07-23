function [rotmat,P]=rotmat_proj_init(nterms,beta,alpha)
%ROTMAT_PROJ_INIT Initialize the rotation matrices via stable projection scheme.
%
%  [ROTMAT,P] = rotmat_proj_init(NTERMS,BETA) initializes 
%  the rotation matrices of degree NTERMS about the y-axis by degree BETA. 
%
%  [ROTMAT,P] = rotmat_proj_init(NTERMS,BETA,ALPHA) initializes 
%  the rotation matrices of degree NTERMS about the z-axis by degree ALPHA 
%  and about the y-axis by degree BETA. 
%
%  After rotation, the expansion pole is moved to location (beta, alpha) 
%  in spherical coordinates (theta, phi).  
%
%  ROTMAT is (NTERMS+1)-by-(2*NTERMS+1)-by-(2*NTERMS+1) complex matrix.
%
%
%  Fast and stable algorithm for initializing the rotation operator about
%  the y-axis determined by angle beta. 
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

if( nargin == 0 ), rotmat_proj_init_test(); return; end;

if( nargin < 3 ), alpha = 0; end;

%%%rotmat = zeros(nterms+1,2*nterms+1,2*nterms+1);
rotmat = zeros((nterms+1)*(2*nterms+1)*(2*nterms+1),1);

ier = 0;
mex_id_ = 'rotmat_projf90(i double[x], i double[x], i int[x], io dcomplex[], i int[x])';
[rotmat] = rotgrid_r2014a(mex_id_, beta, alpha, nterms, rotmat, nterms, 1, 1, 1, 1);

if( alpha == 0 ), rotmat = real(rotmat); end;

if( nargout > 1 ),
    rotmat = reshape(rotmat,nterms+1,2*nterms+1,2*nterms+1);
    P = cell(nterms+1,1);
    for i=1:nterms+1,
       ind = nterms+1-i+1:nterms+1+i-1;
       P{i}=squeeze(rotmat(i,ind,ind));
    end
end

function rotmat_proj_init_test()

nterms = 2
beta = pi/2
a=real(rotmat_proj_init(nterms,beta))

for i=1:nterms+1,
  squeeze(a(i,nterms+1-i+1:nterms+1+i-1,nterms+1-i+1:nterms+1+i-1))
end


