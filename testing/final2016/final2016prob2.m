clear all; clc; format long;
%% a & b)
X0=[3;0;2.1;0];
h=0.01;
N=60/h-1;

[X,t]=RK4(X0,N,h,@fun);
figure
plot(t(1:5/h),X(1,1:5/h),'r'); grid minor;
xt=X(1,:)';
xw=dftmat(xt);
fc=abs(xw);
N=length(xt);
k=[-N/2:1:N/2-1];
wk=((2*pi)/(N*h)).*k;
figure
semilogy(wk,fc,'.m');grid minor;
[fcsorted,I]=sort(fc);
I=I(end-4:end,1);
wkmax=wk(I)
% c)
x=X(1,:);
px=X(2,:);
y=X(3,:);
py=X(4,:);

H=0.5*(px.^2+py.^2)+(x.^2+y-11).^2+(x+y.^2-7).^2;
var=abs(H-H(1));
figure
semilogy(t,var,'.r'); grid minor;
X0=[3 0 2 0]';
DH=jacmod(@fun2,X0,0);
format rat
DH
format long
[veps,vaps]=eig(DH);
vaps=diag(vaps);
vaps
