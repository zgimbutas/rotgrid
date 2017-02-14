%
% Test interpolation from Legendre to Chebychev nodes and back...
%

addpath '../matlab';
addpath '../nufft1d';
addpath '../nufft2d';


n=12*2;

[xc,wc]=chebexps(n);
[xg,wg]=legeexps(n);

a=legematrin(n,n,xc);
cond(a)

b=chebmatrin(n,n,xg);
cond(b)

error=norm(a*b-eye(n,n))

yc=sin(2*pi*xc);
yg=sin(2*pi*xg);

yg' - b*yc'
yc' - a*yg'

if( 1 == 2 ),
theta = 2*pi*((1:n)-0.5)/n
c = exp(-1i*kron(0:n-1,theta'));

[ts]=grule(n/2);
ts=fliplr(ts);
theta=[acos(ts) 2*pi-fliplr(acos(ts))];
d = exp(-1i*kron(0:n-1,theta'));
end
