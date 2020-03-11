function x = qrSolveAlvaro(UR, b);
    [n, m] = size(UR);
    if n > m; M = m; else M = m - 1; end

    for jj = 1:M
        u = [1; UR(jj + 1:n, jj)]; gam = 2 / norm(u)^2;
        b(jj:n) = b(jj:n) - gam * u * (u' * b(jj:n));
    end

    x = BS(triu(UR(1:m, 1:m)), b(1:m));
