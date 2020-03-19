function integral = cuadratura_tanh(a, b, c, n, fun)
%TANH Summary of this function goes here
%   Cuadratura utilitzant el metode tanh. Mes rapid que clenshaw i fejer.

h = c/(n)^(1/2);
sumatori = 0;
for j = [-n:1:n]
    sumatori = sumatori + fun(tanh((j*h)/2))*(cosh((j*h)/2))^(-2);
end

integral = (h/2)*sumatori;

end

