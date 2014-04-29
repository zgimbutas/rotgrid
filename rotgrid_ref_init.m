function fmodes=rotgrid_ref_init(nphi,phi,ntheta,theta,fgrid)
%ROTGRID_REF_INIT Get Fourier modes for a spherical grid.
%
%  Get Fourier modes for a function 
%  defined on a double Chebychev spherical grid.
%
%
%  The output is NPHI-by-NTHETA complex matrix, which contains
%  the values of the Fourier coefficients. Note, that the user data
%  must be specified as NPHI-by-NTHETA complex matrix.
%
%
%  Input parameters:
%
%  nphi - the number of points in lattitude discretization
%  phi - lattitude discretization angles
%  ntheta - the number of points in big circle discretization
%  theta - big circle discretization angles
%  fgrid - user data, NPHI-by-NTHETA complex matrix
%
%  Output parameters:
%
%  fmodes - Fourier modes, NPHI-by-NTHETA complex matrix
%
  
eps=1e-13;


%
% Optionally, use nuFFT in R^2 to compare timings
%
if_use_nufft2d = 0;
if( if_use_nufft2d == 1 ),
%
% Get Fourier coefficients for lattitudes
% Get Fourier coefficients for great circles
iflag=-1;

phi2=repmat(phi,1,ntheta);
theta2=repmat(theta,nphi,1);

fmodes=nufft2d1(nphi*ntheta,phi2,theta2,fgrid,iflag,eps,nphi,ntheta);

return
end


%
% Hybrid method
% Use nuFFT/FFT 1d in theta, and nuFFT/FFT 1d in phi
%

if_use_nufft1d_phi=0;
if( if_use_nufft1d_phi == 1 ),
fgrid_hat=zeros(nphi,ntheta);
%
% Get Fourier coefficients for lattitudes, for arbitrary spaced phi
%
fgrid_hat=zeros(nphi,ntheta);
iflag=-1;
for k=1:ntheta,
  fgrid_hat(:,k)=nufft1d1(nphi,phi,fgrid(:,k),iflag,eps,nphi);
end
else
%
% Get Fourier coefficients for lattitudes, for uniformly spaced phi
%
fgrid_hat=fftshift(fft(fgrid),1)/nphi;
end


if_use_nufft1d_theta=0;
if( if_use_nufft1d_theta == 1 ),
fmodes=zeros(ntheta,nphi);
%
% Get Fourier coefficients for great circles, for arbitrary spaced theta
%
fmodes=zeros(nphi,ntheta);
iflag=-1;
for k=1:nphi,
  fmodes(k,:)=nufft1d1(ntheta,theta,fgrid_hat(k,:),iflag,eps,ntheta);
end
else
%
% Get Fourier coefficients for great circles, for uniformly spaced theta
%
fmodes = fftshift(fft(fgrid_hat,[],2),2)/ntheta;
if( mod(ntheta,2) == 0 ),  
  phase_shift = exp(-1i*(-ntheta/2:(ntheta-1)/2)*.5/ntheta*2*pi);
else
  phase_shift = exp(-1i*(-(ntheta-1)/2:(ntheta-1)/2)*.5/ntheta*2*pi);
end
for k=1:nphi,
  fmodes(k,:) = fmodes(k,:) .* phase_shift;
end
end
