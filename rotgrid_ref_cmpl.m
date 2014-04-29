function grids=rotgrid_ref_cmpl(nphi,phi,ntheta,theta,fmodes,ngrid,xgrid,ygrid,zgrid)
%ROTGRID_REF_CMPL: Fast rotation of a grid by nphi uniformly spaced angles.
%
%  Fast rotation of the user-defined grid (xgrid,ygrid,zgrid).
%  This function constructs a set of grids, obtained by rotating 
%  the original grid by nphi uniformly spaced angles phi.
%
%  The output is NGRID-by-NPHI complex matrix, containing the rotated grids.
%
%  Reference Matlab code...
%
%
%  Input parameters:
%
%  nphi - the number of points in lattitude discretization
%  phi - lattitude discretization angles
%  ntheta - the number of points in big circle discretization
%  theta - big circle discretization angles
%  fmodes - Fourier modes, NPHI-by-NTHETA complex matrix
%  ngrid - the number of points in the user-defined grid
%  xgrid - the x-coordinates of the grid
%  ygrid - the y-coordinates of the grid
%  zgrid - the z-coordinates of the grid
%
%  Output parameters:
%
%  grids - function values at rotated grids, NGRID-by-NPHI complex matrix
%

eps=1e-13;


%
% ...and the angles to be interpolated are
%
rgrid = sqrt(xgrid.^2+ygrid.^2);
zang = atan2(rgrid,zgrid);

%
% ...and the phases to be interpolated are
%
phase = atan2(ygrid,xgrid);


%
% Unwrap phases to interval [-pi,pi], if needed
%
if( 1 == 2 ),
for i=1:nphi,
for j=1:ntheta,
  if( zang(i,j) > pi ), zang(i,j)=zang(i,j)-2*pi; end
end
end
end


%
% Optionally, use nuFFT in R^2 to compare timings
%
if_use_nufft2d = 0;
if( if_use_nufft2d == 1 ),
%
% Call nonuniform FFT in R^2 to get all grids simultaneously
%
xii=reshape(zang,ngrid,1);
pii=reshape(phase,ngrid,1);

grid_thetas=repmat(xii,1,nphi);
grid_phis=repmat(phi,ngrid,1)+repmat(pii,1,nphi);

iflag=+1;
grids=nufft2d2(ngrid*nphi,grid_phis,grid_thetas,iflag,eps,nphi,ntheta,fmodes);
grids=reshape(grids,ngrid,nphi);

return
end


%
% Hybrid method
% Use nuFFT 1d in theta, and nuFFT/FFT 1d in phi
%

%
% Call nonuniform FFT to get all lattitudes
%
xii=reshape(zang,ngrid,1); 
ifsort = 0;
if( ifsort == 1 ), [xii,i]=sort(xii); end
cii=zeros(ngrid,nphi);
%
% Get interpolated values in theta
iflag=+1;
if( 1 == 2 ),
% regular nufft1d2
for k=1:nphi,
  cii(:,k)=nufft1d2(ngrid,xii,iflag,eps,ntheta,fmodes(k,:));
  if( ifsort == 1 ), cii(i,:)=cii(1:ngrid,:); end
end
else
% vectorized nufft1d2 
  cii=nufft1d2v(nphi,ngrid,xii,iflag,eps,ntheta,fmodes.');
  if( ifsort == 1 ), cii(i,:)=cii(1:ngrid,:); end
end


%
% Shift the grids on lattitudes (apply phase in Fourier domain)
%
pii=reshape(phase,ngrid,1);
grids=zeros(ngrid,nphi);

if_use_nufft1d_phi=0;
if( if_use_nufft1d_phi == 1 ),
iflag=+1;
for i=1:ngrid
  grids(i,:) = nufft1d2(nphi,phi+pii(i),iflag,eps,ntheta,cii(i,:));
end
else
%if( mod(nphi,2) == 0 ),  
%for i=1:ngrid
%  cii(i,:) = cii(i,:) .* exp(+1i*((-nphi/2:(nphi-1)/2))*pii(i));
%end
%else
%for i=1:ngrid
%  cii(i,:) = cii(i,:) .* exp(+1i*((-nphi/2:(nphi-1)/2)+.5)*pii(i));
%end
%end
cii=shiftphase(cii,ngrid,nphi,pii);
grids = ifft(ifftshift(cii,2),[],2) *nphi;
end


