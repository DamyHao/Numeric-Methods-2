function [Y]=newtonfun(X)
theta=X(1,1);
T=X(2,1);
N=20000;
h=T/N;


x0=-3.779310253377747;
vx0=17*cos(theta);
y0=-3.283185991286169;
vy0=17*sin(theta);

X0=[x0;vx0;y0;vy0];
[X,t]=RK4(X0,N,h,@funb);
xf=X(1,end);
yf=X(3,end);
xc=-2.805118086952745;  yc=3.131312518250573;
Y=[xf-xc;yf-yc];
end