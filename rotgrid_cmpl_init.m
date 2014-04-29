function fmodes=rotgrid_cmpl_init(nphi,phi,ntheta,theta,fgrid)
fmodes = zeros(nphi,ntheta)+1i*zeros(nphi,ntheta);

mex_id_ = 'rotgridi(i int[], i double[], i int[], i double[], i dcomplex[], io dcomplex[])';
[fmodes] = rotgrid_r2014a(mex_id_, nphi, phi, ntheta, theta, fgrid, fmodes);


