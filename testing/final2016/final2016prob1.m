clear all; clc; format long g;
x0=[-6:1:6];
y0=[-6:1:6];
tol=10-14; itmax=100;crit=[];
for i=1:length(x0)
    for j=1:length(y0)
        X0=[x0(i); y0(j)];
        [extk,res,it]=newtonNLS(X0,tol,itmax,@jacV);
        ext=extk(:,end);
        crit=[crit ext];
    end
end
crit=crit';
C=uniquetol(crit,10^-6,'ByRows',true);
C=C';
figure
plot(C(1,:),C(2,:),'.r');grid minor; axis([-6 6 -6 6]);

sad=[];max=[];min=[];
for n=1:length(C(1,:))
    H=hessV(C(:,n));
   if det(H)<0
       fprintf(' The point \n x= %d \n y=%d\n is a saddle\n\n',C(1,n),C(2,n));
       sad=[sad n];
   else
       if H(1,1)<0 && H(2,2)<0
         fprintf(' The point \n x= %d \n y=%d\n is a maximum\n\n',C(1,n),C(2,n));  
         max=[max n];
       else
         fprintf(' The point \n x= %d \n y=%d\n is a minimum\n\n',C(1,n),C(2,n));  
         min=[min n];
       end
   end
end
xx=[-5:0.1:5];
yy=[-5:0.1:5];
Z=[];
for i=1:length(xx)
    for j=1:length(yy)
        Z(i,j)=V(xx(i),yy(j));
    end
end
Z=Z';

tol=10^-14; itmax=100;
x01L=1;y01L=-1;
x02L=1;y02L=-1.01;

x01S=-2.5;y01S=-2;
x02S=-2.5;y02S=-2.01;

[XK1L,res,it]=newtonNLSmod(x01L,tol,itmax,@level,y01L);
[XK2L,res,it]=newtonNLSmod(x02L,tol,itmax,@level,y02L);
[XK1S,res,it]=newtonNLSmod(x01S,tol,itmax,@level,y01S);
[XK2S,res,it]=newtonNLSmod(x02S,tol,itmax,@level,y02S);

Y1L=[XK1L(:,end);-1];
Y2L=[XK2L(:,end);-1.01];
Y1S=[XK1S(:,end);-2];
Y2S=[XK2S(:,end);-2.01];
YL=[Y1L Y2L];
YS=[Y1S Y2S];
tol=10^-8; itmax=1000;
for i=1:3000
[yl,itconv]=contmod(@level2,YL(:,i),YL(:,i+1),1,tol,itmax);
[ys,itconv]=contmod(@level2,YS(:,i),YS(:,i+1),1,tol,itmax);
YL=[YL yl];
YS=[YS ys];
end

figure
surf(xx,yy,Z);xlabel('x');ylabel('y');zlabel('V(x,y)');title('V(x,y)');hold on;
colormap(hsv)
colorbar
for n=1:length(C(1,:))
    potential=V(C(1,n),C(2,n));
    plot3(C(1,n),C(2,n),potential,'.c','MarkerSize',30);
end
VVl=100*ones(1,length(YL(1,:)));
VVs=100*ones(1,length(YS(1,:)));
plot3(YL(1,:),YL(2,:),VVl,'.m');
plot3(YS(1,:),YS(2,:),VVs,'.m');