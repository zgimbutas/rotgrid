function cii=shiftphase(cii,nii,nphi,pii)
%SHIFTPHASE: Apply a phase shift to each row
%
%   cii=SHIFTPHASE(cii,nii,nphi,pii);
%
%    cii(i,:) = cii(i,:) .* exp(+1i*(-nphi/2:(nphi-1)/2)*pii(i));
%
%  Input parameters:
%
%     cii     Fourier coefficients, 
%             dimensioned as (nii,nphi) array (complex *16)
%     pii     Phase factors for each row (nii) (real *8)     
%                 
%  Output parameters:
%
%     cii     Fourier coefficients (complex *16)
%
%

mex_id_ = 'shiftphase(io dcomplex[], i int[x], i int[x], i double[])';
[cii] = rotgrid_r2014a(mex_id_, cii, nii, nphi, pii, 1, 1);


