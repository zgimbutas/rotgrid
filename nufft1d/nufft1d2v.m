function [cj,ier]=nufft1d2v(nvec,nj,xj,iflag,eps,ms,fk)
%NUFFT1D2V: Nonuniform FFT in R^1 - Type 2, vectorized
%
%  [CJ,IER] = NUFFT1D2V(NVEC,NJ,XJ,IFLAG,EPS,MS,FK);
%
%     cj(j,l) = SUM   fk(k1,l) exp(+/-i k1 xj(j)) 
%                k1  
%                            for j = 1,...,nj
%                            for l = 1,...,nvec
%
%     where -ms/2 <= k1 <= (ms-1)/2
%
%     If (iflag .ge.0) the + sign is used in the exponential.
%     If (iflag .lt.0) the - sign is used in the exponential.
%
%  Input parameters:
%
%     nvec   number of vectors to transform  (integer)
%     nj     number of output values   (integer)
%     xj     location of output values (real *8 array)
%     iflag  determines sign of FFT (see above)
%     eps    precision request  (between 1.0e-15 and 1.0e-1)
%     ms     number of Fourier modes given  [ -ms/2: (ms-1)/2 ]
%     fk     Fourier coefficient values (complex *16 array)
%
%  Output parameters:
%
%     cj     output values (complex *16 array)
%     ier    error return code   
%            ier = 0  => normal execution.
%            ier = 1  => precision eps requested is out of range.
%
%

cj=zeros(nj,nvec)+1i*zeros(nj,nvec);
ier=0;

mex_id_ = 'nufft1d2vf90(i int[x], i int[x], i double[], io dcomplex[], i int[x], i double[x], i int[x], i dcomplex[], io int[x])';
[cj, ier] = nufft1d(mex_id_, nvec, nj, xj, cj, iflag, eps, ms, fk, ier, 1, 1, 1, 1, 1, 1);


