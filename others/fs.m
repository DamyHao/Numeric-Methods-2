function x = fs(L, b)
    % ALVARO
    % Resol un sistema lower trinagular Lx = b
    % b ha de ser una columna
    % Output x: solution of Lx = b
    x = 0 * b; n = length(b); x(1) = b(1) / L(1, 1);

    % Formula:x(i) 1/l(i,i)*{bi - sumatori[l(i,j)*x(j)]} Sumatori de j = 1 a i-1
    for ii = 2:n
        x(ii) = (b(ii) - L(ii, 1:ii - 1) * x(1:ii - 1)) / L(ii, ii);
    end

end
