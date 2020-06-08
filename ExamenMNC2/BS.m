function x = BS(U, b)
    % Resol un sistema upper trinagular Ux = b
    % b ha de ser una columna
    % Output x: solution of Lx = b
    % Formula:x(i) 1/u(i,i)*{b(i) - sumatori[u(i,j)*x(j)]} Sumatori de j = i+1 a n
    [dim, n] = size(U);
    x = zeros(n, 1);

    if dim ~= n; error('NECESSITA UNA MATRIU QUADRADA. Revisar els arguments');

    else
        for i = n:-1:1
            sum = 0;

            for j = i + 1:n
                sum = sum + U(i, j) * x(j);
            end

            x(i) = (1 / U(i, i)) * (b(i) - sum);
        end

    end

end
