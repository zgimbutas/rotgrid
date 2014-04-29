function grids=rotgrid_real_opt_vec(nphi,phi,ntheta,theta,fmodes,ngrid,nlat,xgrid,ygrid,zgrid)
%ROTGRID_REAL_OPT_VEC: Fast rotation of a grid by nphi uniformly spaced angles.
%
%  Fast rotation of the user-defined grid (xgrid,ygrid,zgrid).
%  This function constructs a set of grids, obtained by rotating 
%  the original grid by nphi uniformly spaced angles phi.
%
%  The output is NGRID-by-NPHI-by-NLAT real matrix, 
%  containing the rotated grids.
%
%  Optimized for symmetric spherical grid, nphi is even, ngrid = nphi*ntheta
%  Note, this routine does not work for arbitrary grids, use rotgrid_real.
%
%  Input parameters:
%
%  nphi - the number of points in lattitude discretization
%  phi - lattitude discretization angles
%  ntheta - the number of points in big circle discretization
%  theta - big circle discretization angles
%  fmodes - Fourier modes, NPHI-by-NTHETA complex matrix
%  ngrid - the number of points in the user-defined grid, see note above
%  nlat - the number latitudes/grids
%  xgrid - the x-coordinates of the grid, NGRID-by-NLAT
%  ygrid - the y-coordinates of the grid, NGRID-by-NLAT
%  zgrid - the z-coordinates of the grid, NGRID-by-NLAT
%
%  Output parameters:
%
%  grids - function values at rotated grids, NGRID-by-NPHI-by-NLAT real matrix
%

grids = zeros(ngrid,nphi,nlat);

mex_id_ = 'rotgrid_real_opt_vec(i int[], i double[], i int[], i double[], i dcomplex[], i int[], i int[], i double[], i double[], i double[], io double[])';
[grids] = rotgrid_r2014a(mex_id_, nphi, phi, ntheta, theta, fmodes, ngrid, nlat, xgrid, ygrid, zgrid, grids);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

