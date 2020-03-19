func = @(x)((Vyo/Vxo).*x)-((1/2).*g.*(x.^2)/(Vxo.^2)) - (a*(x.^2).*exp((-b).*x));
[xk,it,p] = newton(70, 80, 10^-3, 200, func);