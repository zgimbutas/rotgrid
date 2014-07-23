nj = 1000;
ms = 2100;

nvec = 10;

xj = sort((rand(nj,1)*2-1)*pi);
cj = randn(nj,nvec)+1i*randn(nj,nvec);

eps=1e-12;


iflag = +1;

tic
fk = zeros(ms,nvec);
for i=1:nvec
fk(:,i) = dirft1d1(nj,xj,cj(:,i),iflag,ms);
end
toc

tic
fk1 = nufft1d1v(nvec,nj,xj,cj,iflag,eps,ms);
toc

abs_error=norm(fk-fk1,2)
rel_error=norm(fk-fk1,2)/norm(fk,2)


iflag = -1;

tic
fk = zeros(ms,nvec);
for i=1:nvec
fk(:,i) = dirft1d1(nj,xj,cj(:,i),iflag,ms);
end
toc

tic
fk1 = nufft1d1v(nvec,nj,xj,cj,iflag,eps,ms);
toc

abs_error=norm(fk-fk1,2)
rel_error=norm(fk-fk1,2)/norm(fk,2)
