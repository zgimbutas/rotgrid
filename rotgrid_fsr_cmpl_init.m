function [rotmat]=rotgrid_fsr_cmpl_init(nterms,nbeta,beta)
%ROTGRID_FSR_CMPL_INIT Initialize rotgrid_fsr_cmpl.
%
%  Fast, FFT-based algorithm for rotating complex spherical harmonic grids
%  into pole locations (beta_j, alpha_k), 
%  where alpha_k = 2*pi * k/nrot, k=0..nrot-1, beta_j=1..nbeta
%
%  Input parameters:
%
%  nterms - the number of terms in spherical harmonics expansion
%  nbeta - the number of points in pole location meridian discretization
%  beta - angles for new pole locations beta_j, j=1..nbeta
%
%  Output parameters:
%
%  rotmat - The rotation operators for directions (beta_j,0). 

rotmat=zeros((nterms+1)*(2*nterms+1)*(2*nterms+1),nbeta);

mex_id_ = 'rotgrid_fsr_cmpl_init(i int[], i int[], i double[], io double[])';
[rotmat] = rotgrid_r2014a(mex_id_, nterms, nbeta, beta, rotmat);


