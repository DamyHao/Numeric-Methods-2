function [Y] = funcioAlfa(x, alfa)
g = 9.81;
a = 0.15;
b = 0.04;

Vo = 37;

Vxo = Vo*cosd(alfa);
Vyo = Vo*sind(alfa);

Ypar = ((Vyo/Vxo).*x)-((1/2).*g.*(x.^2)/(Vxo.^2));
Yper = a*(x.^2).*exp((-b).*x);
Y = Ypar - Yper;