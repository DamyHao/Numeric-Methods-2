function [M1, M2, M3, Lhat, x] = crearMatriusODE(n, C, p, q, r, a, b)
    % Input: n: nombre de nodes Chebyshev,
    % C: Matriu 2*3 amb les Robin conditions
    % Funcions p, q, r:  En cas que sigui una funcio constant, cal que apareixi la x en el function handler

    % Fem la diferenciacio Chebyshev, que ens retorna també els nodes:
    [D, x] = chebdiff(n, a, b);
    P = diag(p(x));
    Q = diag(q(x));
    R = diag(r(x));
    L = P * D^2 + Q * D + R;

    %Cas en que nomes tenim una condicio de Robin corresponent a x=1. Haurem de
    %treure el primer element ja que els nodes chebyshev estan al reves.

    Lhat = L(2:end, 2:end); %treiem primera fila i primera columna de L
    % En el cas normal es treu tambe ultima fila i ultima
    M1 = -D(1, 2:end); %Primera fila de D, sense agafar element de la primera columna

    M2 = [C(2, 1) + C(2, 2) * D(1, 1)];

    M3 = [L(2:end, 1)];

end
