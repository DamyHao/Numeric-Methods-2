function I_fej= fejer(a, b, N, fun)
%   Cuadratura de Fejer, que es com la de Clenshaw pero oberta.
%   La funcio entregada fun ha de ser capaç de tractar vectors si se li
%   donen com arguments.

if(1 == mod(N, 2))
    error('"N" HA DE SER PARELL!!');
end

k = 1:N;
Wj = zeros(1,N);
thk = ((2.*k-1).*pi)/(2.*N);
posicio = 0;

for j = thk
    posicio = posicio + 1;
    sumatori=0;
    for i = 1:N/2
        sumatori = sumatori + (cos(2*i*j)/(4*(i^2)-1));
    end
    
    Wj(posicio) = (2/N)*(1-2*sumatori);
end

xk = cos(thk);
phi = b + (b-a).*(1/2).*(xk-1);
f_phis = fun(phi);
I_fej = ((b-a)/2).*(Wj*f_phis');
end

