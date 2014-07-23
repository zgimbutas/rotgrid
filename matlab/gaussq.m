function [t,w] = gaussq(kind,n,alpha,beta,kpts,endpts)

b = zeros(n,1);
t = zeros(n,1);
w = zeros(n,1);

mex_id_ = 'gaussq(i int[x], i int[x], i double[x], i double[x], i int[x], i double[x], io double[], io double[], io double[])';
[b, t, w] = rotgrid_r2014a(mex_id_, kind, n, alpha, beta, kpts, endpts, b, t, w, 1, 1, 1, 1, 1, 2);


