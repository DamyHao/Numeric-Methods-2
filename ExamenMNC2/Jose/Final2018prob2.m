%% a)
clear all;clc;format long;
x0A=[-4;-4];x0B=[-4; 0];x0C=[-4; 4];
tol=10^-15;itmax=100;


[XA,resA,itA]=newtonNLS(x0A,tol,itmax,@funa);
[XB,resA,itA]=newtonNLS(x0B,tol,itmax,@funa);
[XC,resA,itA]=newtonNLS(x0C,tol,itmax,@funa);
XA(:,end)
XB(:,end)
XC(:,end)
%% b)

h=0.01; tf=60; N=tf/h+1;
f0=[-2.8;0;3.1;0];
[f,t]=RK4(f0,N,h,@funb);
figure
plot(t(1:501),f(1,1:501));grid minor;xlabel('t');ylabel('x(t)');title('x(t) for t\in[0,5]');
figure
plot(t(1,:),f(1,:));grid minor;xlabel('t');ylabel('x(t)');title('x(t) for t\in[0,60]');
% c)
xt=f(1,:)';
NN=N+1;
fsf=dftmat(xt);
k=[-NN/2:1:NN/2-1];
wk=((2*pi)/(NN*h))*k;
figure
semilogy(wk,abs(fsf),'-m');grid minor; xlabel('\omega_k');ylabel('|fourier coeficients|');title('c fourier analysis of x(t)');

%% d)
% from a) we know (x,y) coordinates from C are xc=-2.805118086952745  yc=3.131312518250573
xc=-2.805118086952745;  yc=3.131312518250573;
X0=[xc 0 yc 0]'; %(the momentums are 0 since C is an equilibrium point);
Df=jacmod(@funb2,X0,t);
eig(Df)


%% e)
clear all;clc;
tol=10^-8;itmax=100;

x01=[0.7;0.8];x02=[1.13;0.83];x03=[1.6;0.55];
[Xk1,res,it]=newtonNLS(x01,tol,itmax,@newtonfun);
[Xk2,res,it]=newtonNLS(x02,tol,itmax,@newtonfun);
[Xk3,res,it]=newtonNLS(x03,tol,itmax,@newtonfun);

theta1=Xk1(1,end);theta2=Xk2(1,end);theta3=Xk3(1,end);
T1=Xk1(2,end);T2=Xk2(2,end);T3=Xk3(2,end);
N=20000;
h1=T1/N;h2=T2/N;h3=T3/N;
X01=[-3.779310253377747;17*cos(theta1); -3.283185991286169;17*sin(theta1)];
X02=[-3.779310253377747;17*cos(theta2); -3.283185991286169;17*sin(theta2)];
X03=[-3.779310253377747;17*cos(theta3); -3.283185991286169;17*sin(theta3)];

[X1,t]=RK4(X01,N,h1,@funb);
[X2,t]=RK4(X02,N,h2,@funb);
[X3,t]=RK4(X03,N,h3,@funb);
figure
plot(X1(1,:),X1(3,:),'r');grid minor;axis equal; hold on;
plot(X2(1,:),X2(3,:),'g');
plot(X3(1,:),X3(3,:),'b');
plot(-3.779310253377747,-3.283185991286169,'ok');
plot(-2.805118086952745,3.131312518250573,'ob');







