%
% Rotation of spherical harmonics expansions.
% -----------------------------------------------------
%
% rotviaproj_cmpl - Rotate the complex spherical harmonics expansion 
%              about the y-axis via pseudo-spectral projection.
% rotviaproj_real - Rotate the real spherical harmonics expansion 
%              about the y-axis via pseudo-spectral projection.
%
% rotviarecur_cmpl - Rotate the complex spherical harmonics expansion
%              about the y-axis using recurrence.
% rotviarecur_real - Rotate the real spherical harmonics expansion
%              about the y-axis using recurrence.
%
%
% Rotation of spherical harmonics expansions.
% -----------------------------------------------------
%
% rotmat_proj_init - Initialize the rotation matrices
%              about the y-axis via pseudo-spectral projection.
%
% rot1lat_fsr_cmpl - Rotate the complex spherical harmonics expansion
%              into a collection of new pole locations (beta, alpha_k)
%              using Wigner rotation matrices + FFT.
%
% rot1lat_fsr_real - Rotate the real spherical harmonics expansion
%              into a collection of new pole locations (beta, alpha_k)
%              using Wigner rotation matrices + FFT.
%
% rot1lat_cmpl_proj - Rotate the complex spherical harmonics expansion 
%              into a collection of new pole locations (beta, alpha_k)
%              via pseudo-spectral projection.
%
%
% Rotation of spherical harmonics grids.
% Direct rotation + spherical transform
% -----------------------------------------------------
%
% rotgrid_dsr_cmpl - Rotate the complex spherical grid
%              into a collection of new pole locations (beta_j, alpha_k)
%              using direct rotation + spherical transform
% rotgrid_dsr_real - Rotate the real spherical grid
%              into a collection of new pole locations (beta_j, alpha_k)
%              using direct rotation + spherical transform
%
%
% Rotation of spherical harmonics grids.
% Wigner rotation matrices + FFT + spherical transform.
% -----------------------------------------------------
%
% rotgrid_fsr_cmpl - Rotate the complex spherical grid
%              into a collection of new pole locations (beta_j, alpha_k)
%              using Wigner rotation matrices + FFT + spherical transform
% rotgrid_fsr_real - Rotate the real spherical grid
%              into a collection of new pole locations (beta_j, alpha_k)
%              using Wigner rotation matrices + FFT + spherical transform
%
%
% Rotation of spherical harmonics grids.
% Hybrid nuFFT interpolation scheme.
% -----------------------------------------------------
%
% rotgrid_init_dcheb - Get Fourier modes for a double spherical Chebychev grid.
%
% rotgrid_cmpl - Rotate the complex spherical grid
%           into a collection of new pole locations (beta_j, alpha_k)
%           using hybrid nuFFT interpolation scheme.
% rotgrid_real - Rotate the real spherical grid
%              into a collection of new pole locations (beta_j, alpha_k)
%              using hybrid nuFFT interpolation scheme.
%
% rotgrid_cmpl_opt - Rotate the complex spherical grid
%           into a collection of new pole locations (beta_j, alpha_k)
%           using hybrid nuFFT interpolation scheme.
%           Assumes a symmetric spherical grid, nphi is even.
% rotgrid_real_opt - Rotate the real spherical grid
%              into a collection of new pole locations (beta_j, alpha_k)
%              using hybrid nuFFT interpolation scheme.
%              Assumes a symmetric spherical grid, nphi is even.
%
% rotgrid_cmpl_opt_vec - Rotate the complex spherical grid
%           into a collection of new pole locations (beta_j, alpha_k)
%           using hybrid nuFFT interpolation scheme.
%           Assumes a symmetric spherical grid, nphi is even.
%           Vectorized for multiple latitudes/grids.
% rotgrid_real_opt_vec - Rotate the real spherical grid
%              into a collection of new pole locations (beta_j, alpha_k)
%              using hybrid nuFFT interpolation scheme.
%              Assumes a symmetric spherical grid, nphi is even.
%              Vectorized for multiple latitudes/grids.
%
%
% Rotation of spherical harmonics grids.
% Hybrid nuFFT interpolation scheme.  Reference Matlab code.
% -----------------------------------------------------
%
% rotgrid_ref_init - Get Fourier modes for a double spherical Chebychev grid.
%
% rotgrid_ref_cmpl - Rotate the complex spherical grid
%           into a collection of new pole locations (beta_j, alpha_k)
%           using hybrid nuFFT interpolation scheme.
%
%
% Utility functions.
% ------------------
%
% fftnext235 - Returns the next multiple of 2, 3, and 5.
%
%
% Spherical transforms.
% ---------------------
%
% sphtrans_cmpl - Complex spherical transform on a spherical grid.
%        Inverse transform, convert spherical harmonic expansion to grid values.
% sphtrans_real - Real spherical transform on a spherical grid.
%        Inverse transform, convert spherical harmonic expansion to grid values.
%
% sphtrans_fwd_cmpl - Complex spherical transform on a spherical grid.
%        Forward transform, convert grid values to spherical harmonic expansion.
% sphtrans_fwd_real - Real spherical transform on a spherical grid.
%        Forward transform, convert grid values to spherical harmonic expansion.
%
% sphtrans_cmpl_lege_init - Initialize complex spherical transforms 
%                           on a spherical Legendre grid.
% sphtrans_cmpl_cheb_init - Initialize complex spherical transforms
%                           on a spherical Chebychev grid.
% sphtrans_real_lege_init - Initialize real spherical transforms
%                           on a spherical Legendre grid.
% sphtrans_real_cheb_init - Initialize real spherical transforms
%                           on a spherical Chebychev grid.
%
%
% Orthogonal polynomials (interpolation and integration).
% -------------------------------------------------------
%
% legeexps - Gauss-Legendre quadrature on the interval [-1,1].
% chebexps - Chebychev quadrature on the interval [-1,1].
%
% legematrin - Interpolation matrix for Legendre nodes.
% chebmatrin - Interpolation matrix for Chebychev nodes.
%
% legecheb_double - Interpolation matrix for double Legendre grids.
%
%
% Spherical grids.
% ----------------
%
% xyz_grid - Construct a spherical grid 
%            and rotate it by angle beta around y-axis.
% xyz_grid_a - Construct a spherical grid 
%            rotate it by angle beta around y-axia, 
%            then, rotate it by angle alpha around z-axis.
%
% init_grid_lege_single - Parameters for single spherical Legendre grid.
% init_grid_lege_double - Parameters for double spherical Legendre grid.
% init_grid_cheb_single - Parameters for single spherical Chebychev grid.
% init_grid_cheb_double - Parameters for double spherical Chebychev grid.
%
% Testing routines.
% -----------------
%
% test_rotgrid_dsr* - test direct grid rotation scheme
% test_rotgrid_fsr* - test FFT-accelerated grid rotation scheme
% test_rotgrid_hnufft* - test hybrid nuFFT grid rotation scheme
%
