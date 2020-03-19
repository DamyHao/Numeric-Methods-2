function [y] = diferencia(x, Vo, A);

g = 9.81; a = 0.15; b = 0.04;

Vxo = Vo*cosd(A);
Vyo = Vo*sind(A);

Ypar = ((Vyo/Vxo).*x)-((1/2).*g.*(x.^2)/(Vxo.^2));
Yper = a*(x.^2).*exp((-b).*x);

y = Ypar - Yper;
end
