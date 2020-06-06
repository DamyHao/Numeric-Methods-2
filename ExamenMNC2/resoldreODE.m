function [f, x] = resoldreODE(n, C, p, q, r, g, a, b)
    % Les funcions entregades han de ser les corresponents al domini -1, 1
    % Totes les funcions han ser aptes per vectors
    [M1, M2, M3, Lhat, x] = crearMatriusODE(n, C, p, q, r,a,b);
    mCoef = [C(2, 2), 0;
        0, C(1, 2)];
    esquerra = Lhat + M3 * inv(M2) * mCoef * M1;
    dreta = g(x(2:end - 1)) - M3 * inv(M2) * [C(2, 3); C(1, 3)];

    f = esquerra \ dreta;
end
