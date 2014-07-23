function f=funuser_real(x,y,z)
%
%  User-defined function
%

%%%f = x.^2+x+sin(y+x)+z + 1i*(x+y.^2+sin(x+y+z.^2));
f = x.^2+x+2*y+z+1 + 1i*(x+2*z+y+y.^2);
%%%f = x.^2;
%%%f = x.^2+x+z + y.^2;

%f = x;
%f = y;
%f = z;

f = real(f);

