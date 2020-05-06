function x = BS(A, b)

[dim, n] = size(A);
x = zeros(n, 1);

if dim ~= n; error('NECESSITA UNA MATRIU QUADRADA. Revisar els arguments');
    
else
    
    for i = n:-1:1
        sum = 0;
        
        for j = i + 1:n
            sum = sum + A(i, j) * x(j);
        end
        x(i) = (1/A(i,i))*(b(i)-sum);
    end
    
end

end
