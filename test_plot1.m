%  Testing routines for fast gridding algorithm
%
%  Fast rotation of the user-defined grid (xrot,yrot,zrot).
%  This function constructs a set of grids, obtained by rotating 
%  the original grid by nphi uniformly spaced angles phi.


addpath '../nufft1d';
addpath '../nufft2d';


%
% Construct the initial spherical grid, 
% it will be used to get Fourier coefficients for a sphere.
%
% NOTE: we are using great circles theta = [-pi..pi] 
% (and not the meridians theta = [0..pi]), in order 
% to simplify the code. 
%

p = 11;
ntheta0=p+1;

ntheta=ntheta0*2;
nphi=ntheta0*2;

nphi=fftnext235(nphi);

[phi,theta,xs,ys,zs]=init_grid_cheb_double(nphi,ntheta);

%
% Get the values of the user-defined function 
% on the initial spherical grid (theta in [-pi..pi], see note above)
%
beta = 0;

[xgrid,ygrid,zgrid]=xyz_grid(beta,nphi,xs,ys,ntheta,zs,theta);
fgrid = funuser_real(xgrid,ygrid,zgrid);

figure(1)
plot3(xgrid(:),ygrid(:),zgrid(:),'*')
xlabel('x');
ylabel('y');
zlabel('z');
axis equal
axis off

%
% Construct the rotated grid. Test an arbitrary ntheta grid.
%
% The algorithm does not assume that the rotated points are
% structured in any way, the user can specify arbitrary number
% and locations of points. In this test, we check a regular
% spherical grid.
%

nphi1=nphi
ntheta1=ntheta/2

[phi1,theta1,xs1,ys1,zs1]=init_grid_single(nphi1,ntheta1);

beta = acos(zs1(3));

[xrot,yrot,zrot]=xyz_grid(beta,nphi1,xs1,ys1,ntheta1,zs1,theta1);
frot = funuser_real(xrot,yrot,zrot);

figure(2)
plot3(xrot(:),yrot(:),zrot(:),'*')
xlabel('x');
ylabel('y');
zlabel('z');
axis equal
axis off


nrot=size(xrot,1)*size(xrot,2)

ngrid= nrot*ntheta1

xgrid= zeros(nphi1,ntheta1,ntheta1);
ygrid= zeros(nphi1,ntheta1,ntheta1);
zgrid= zeros(nphi1,ntheta1,ntheta1);
for i=1:ntheta1
[xrot1,yrot1,zrot1]=xyz_grid(zs1(i),nphi1,xs1,ys1,ntheta1,zs1,theta1);
xgrid(:,:,i)=xrot1(:,:);
ygrid(:,:,i)=yrot1(:,:);
zgrid(:,:,i)=zrot1(:,:);
end


nrot = nphi1
xgrid= zeros(nphi1,ntheta1,nrot);
ygrid= zeros(nphi1,ntheta1,nrot);
zgrid= zeros(nphi1,ntheta1,nrot);
for i=1:nrot
alpha=(i-1)*2*pi/nrot;
beta=acos(zs1(3));
[xrot1,yrot1,zrot1]=xyz_grid_a(beta,nphi1,xs1,ys1,ntheta1,zs1,theta1,alpha);
xgrid(:,:,i)=xrot1(:,:);
ygrid(:,:,i)=yrot1(:,:);
zgrid(:,:,i)=zrot1(:,:);
end


figure(3)
plot3(xgrid(:),ygrid(:),zgrid(:),'.')
xlabel('x');
ylabel('y');
zlabel('z');
axis equal
axis off

figure(1)
print -dpng 'p1.png' -r300
figure(2)
print -dpng 'p2.png' -r300
figure(3)
print -dpng 'p3.png' -r300

if( 1 == 2 ),

for i=1:nphi1
figure(10+i)
x1=xgrid(:,:,i);
y1=ygrid(:,:,i);
z1=zgrid(:,:,i);
plot3(x1(:),y1(:),z1(:),'.')
xlabel('x');
ylabel('y');
zlabel('z');
axis equal
axis off
filename = ['pf',num2str(i),'.png']
print('-dpng',['pf',num2str(i),'.png'],'-r300');
end

end


phase = atan2(ygrid(:,:,1),xgrid(:,:,1));
cosa = cos(phase); sina = sin(phase);

for k=1:nrot

for i=1:nphi1
for j=1:ntheta1
    x1 = xgrid(i,j,k);
    y1 = ygrid(i,j,k);
    z1 = zgrid(i,j,k);
    xgrid1(i,j,k) = x1*cosa(i,j) + y1*sina(i,j);
    ygrid1(i,j,k) = -x1*sina(i,j) + y1*cosa(i,j);
    zgrid1(i,j,k) = z1;
end
end

end


figure(4)
plot3(xgrid1(:),ygrid1(:),zgrid1(:),'.')
xlabel('x');
ylabel('y');
zlabel('z');
axis equal
axis off

figure(4)
print -dpng 'p4.png' -r300
