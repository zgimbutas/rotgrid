function [pol,der]=legepol(x,n)

pol=0;
der=0;
mex_id_ = 'legepol(i double[x], i int[x], i double[x], i double[x])';
rotgrid_r2014a(mex_id_, x, n, pol, der, 1, 1, 1, 1);



