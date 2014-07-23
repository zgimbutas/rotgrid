function pols=legepols(x,n)

n1=n+1;
pols=zeros(n1,1);

mex_id_ = 'legepols(i double[x], i int[x], i double[x])';
rotgrid_r2014a(mex_id_, x, n, pols, 1, 1, n1);



