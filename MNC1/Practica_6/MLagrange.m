function P = MLagrange(z, x);
n = length(x);
m = length (z);
den = [];
num = [];
P = [;];


for i = 1:1:n
    % ens genera el vector fila de la matriu dels polinomis cardinals de
    % lagrange.
       den = x(i)-x;
       den(i) = 1; %así, al multiplicar se queda igual, y así se lo salta.
       for k=1:1:m
           num = z(k)-x;
           num(i)=1;
           P(k,i) = prod(num)/prod(den);
       end
       
end

       