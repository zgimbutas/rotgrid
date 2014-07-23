nj = 1000;
nk = 2100;

nvec = 10;

xj = sort((rand(nj,1)*2-1)*pi);
cj = randn(nj,nvec)+1i*randn(nj,nvec);
sk = sort((rand(nk)*2-1)*pi);

eps=1e-12;


iflag = +1;

tic
fk = zeros(nk,nvec);
for i=1:nvec
fk(:,i) = dirft1d3(nj,xj,cj(:,i),iflag,nk,sk);
end
toc

tic
fk1 = nufft1d3v(nvec,nj,xj,cj,iflag,eps,nk,sk);
toc

abs_error=norm(fk-fk1,2)
rel_error=norm(fk-fk1,2)/norm(fk,2)


iflag = -1;

tic
fk = zeros(nk,nvec);
for i=1:nvec
fk(:,i) = dirft1d3(nj,xj,cj(:,i),iflag,nk,sk);
end
toc

tic
fk1 = nufft1d3v(nvec,nj,xj,cj,iflag,eps,nk,sk);
toc

abs_error=norm(fk-fk1,2)
rel_error=norm(fk-fk1,2)/norm(fk,2)
