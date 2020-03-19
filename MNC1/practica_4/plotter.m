clear all;
g = 9.81;
a = 0.15;
b = 0.04;

x = [0:1/100:200];
Vo = 37;
A = 67.5;

Vxo = Vo*cosd(A);
Vyo = Vo*sind(A);

Ypar = ((Vyo/Vxo).*x)-((1/2).*g.*(x.^2)/(Vxo.^2));
Yper = a*(x.^2).*exp((-b).*x);
D = Ypar - Yper;
%fun = @(x)((Vyo/Vxo).*x)-((1/2).*g.*(x.^2)/(Vxo.^2)) -
%(a*(x.^2).*exp((-b).*x))

plot(x,Ypar,'r');
hold on;
plot(x,Yper,'b');
hold on;
plot(x,D,'g');
axis([0 200 0 60]);