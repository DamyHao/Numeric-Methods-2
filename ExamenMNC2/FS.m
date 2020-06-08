function x = FS(L, b)
% Resol un sistema lower trinagular Lx = b

% Formula:x(i) 1/l(i,i)*{bi - sumatori[l(i,j)*x(j)]} Sumatori de j = 1 a i-1
[dim, n] = size(L);
x = zeros(n, 1);

if dim~= n; error('NECESSITA UNA MATRIU QUADRADA. Revisar els arguments');
    
else
    
    for ii = 1:n
        sum = 0;
        
        for jj = 1:ii-1
            sum = sum + L(ii, jj) * x(jj); % Utilitzem tots els coeficients que ja sabem per calcular el que cal restar a la b
        end
        
        x(ii) = (1/L(ii,ii))*(b(ii)-sum); % Una vegada despejat, sabem la nova incognita.
    end
end
end
