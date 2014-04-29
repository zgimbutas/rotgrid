function grids=rotgrid_cmpl_opt(nphi,phi,ntheta,theta,fmodes,ngrid,xgrid,ygrid,zgrid)
%ROTGRID_CMPL_OPT: Fast rotation of a grid by nphi uniformly spaced angles.
%
%  Fast rotation of the user-defined grid (xgrid,ygrid,zgrid).
%  This function constructs a set of grids, obtained by rotating 
%  the original grid by nphi uniformly spaced angles phi.
%
%  The output is NGRID-by-NPHI complex matrix, containing the rotated grids.
%
%  Optimized for symmetric spherical grid, nphi is even, ngrid = nphi*ntheta
%  Note, this routine does not work for arbitrary grids, use rotgrid.
%
%  Input parameters:
%
%  nphi - the number of points in lattitude discretization
%  phi - lattitude discretization angles
%  ntheta - the number of points in big circle discretization
%  theta - big circle discretization angles
%  fmodes - Fourier modes, NPHI-by-NTHETA complex matrix
%  ngrid - the number of points in the user-defined grid, see note above
%  xgrid - the x-coordinates of the grid
%  ygrid - the y-coordinates of the grid
%  zgrid - the z-coordinates of the grid
%
%  Output parameters:
%
%  grids - function values at rotated grids, NGRID-by-NPHI complex matrix
%

grids = zeros(ngrid,nphi)+1i*zeros(ngrid,nphi);

mex_id_ = 'rotgrid_opt(i int[], i double[], i int[], i double[], i dcomplex[], i int[], i double[], i double[], i double[], io dcomplex[])';
[grids] = rotgrid_r2014a(mex_id_, nphi, phi, ntheta, theta, fmodes, ngrid, xgrid, ygrid, zgrid, grids);


