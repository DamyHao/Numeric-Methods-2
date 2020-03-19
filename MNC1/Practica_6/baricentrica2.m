function B = baricentrica2(z, xj, funxj)
% Funcio feta per nodes chebyshev.!!
n = length(xj);
m = length(z);
den = zeros(length(z),1);
num = zeros(length(z),1);

% Aqui calculem la funció baricentrica en tots els punts del vector z
% utilitzant una interpolacio dels punts xj.
for j = [1:1:n]
    if j == 1 || j == n
        numj = ((-1)^j)*funxj(j)*(z-xj(j)).^(-1).*(1/2);
        num = num + numj';
        denj = ((-1)^j)*(z-xj(j)).^(-1).*(1/2);
        den = den + denj';
    else
        num = num + (((-1)^j)*funxj(j)*(z-xj(j)).^(-1))';
        den = den + (((-1)^j)*(z-xj(j)).^(-1))';
    end
end

B = num./den;

%Aprofita el producte de matrius per fer el sumatori.