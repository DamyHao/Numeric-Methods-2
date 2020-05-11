function x = FS(A, b)

[dim, n] = size(A);
x = zeros(n, 1);

if dim~= n; error('NECESSITA UNA MATRIU QUADRADA. Revisar els arguments');
    
else
    
    for i = 1:n
        sum = 0;
        
        for j = 1:i-1
            sum = sum + A(i, j) * x(j); % Utilitzem tots els coeficients que ja sabem per calcular el que cal restar a la b
        end
        
        x(i) = (1/A(i,i))*(b(i)-sum); % Una vegada despejat, sabem la nova incognita.
    end
end

end
