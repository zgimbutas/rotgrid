function [pol,der]=chebpol(x,n)

pol=0;
der=0;
mex_id_ = 'chebpol(i double[x], i int[x], io double[x], io double[x])';
[pol, der] = rotgrid_r2014a(mex_id_, x, n, pol, der, 1, 1, 1, 1);



