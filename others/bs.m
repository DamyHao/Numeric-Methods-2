function x = bs(U, b)
    % ALVARO
    % Resol un sistema upper trinagular Ux = b
    % b ha de ser una columna
    % Output x: solution of Lx = b
    % Formula:x(i) 1/u(i,i)*{b(i) - sumatori[u(i,j)*x(j)]} Sumatori de j = i+1 a n
    x = 0 * b; n = length(b); x(n) = b(n) / U(n, n);

    for ii = n - 1:-1:1
        x(ii) = (b(ii) - U(ii, ii + 1:n) * x(ii + 1:n)) / U(ii, ii);
    end

end
