nj = 1000;
ms = 2100;

nvec = 10;

xj = sort((rand(nj,1)*2-1)*pi);
fk = randn(ms,nvec)+1i*randn(ms,nvec);

eps = 1e-12;


iflag = +1;

tic
cj = zeros(nj,nvec);
for i=1:nvec
cj(:,i) = dirft1d2(nj,xj,iflag,ms,fk(:,i));
end
toc

tic
cj1 = nufft1d2v(nvec,nj,xj,iflag,eps,ms,fk);
toc
 
abs_error=norm(cj-cj1,2)
rel_error=norm(cj-cj1,2)/norm(cj,2)


iflag = -1;

tic
cj = zeros(nj,nvec);
for i=1:nvec
cj(:,i) = dirft1d2(nj,xj,iflag,ms,fk(:,i));
end
toc

tic
cj1 = nufft1d2v(nvec,nj,xj,iflag,eps,ms,fk);
toc

abs_error=norm(cj-cj1,2)
rel_error=norm(cj-cj1,2)/norm(cj,2)

