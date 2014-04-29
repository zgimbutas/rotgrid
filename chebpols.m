function pols=chebpols(x,n)

n1=n+1;
pols=zeros(n1,1);

mex_id_ = 'chebpols(i double[x], i int[x], io double[x])';
[pols] = rotgrid_r2014a(mex_id_, x, n, pols, 1, 1, n1);



