function [f, x] = resoldreODE(n, C, p, q, r, g)
    % Resol una equacio diferencial de tipus boundary en el interval a,b.
    % ATENCIO: fiquem les C, p, q, r, al domini de (a,b) (en realitat falta per comprovar si va aqui pero 99% segur que si jeje)
    % INPUT:
    %   n: nombre de nodes
    %   C: matriu dels coeficients C
    %       {C_11*f(a) + C_12*f(a) = C_13*f(a)
    %       {C_21*f(b) + C_22*f(b) = C_23*f(b)
    %   p, q, r, g:funcions del boundary problem p(x)f'' + q(x)f' + r(x)f = g(x)*f
    %           En cas que sigui una funcio constant, cal que apareixi la x en
    %           el function handler. Les funcions entregades han de ser les
    %           corresponents al domini a, b.
    %           Totes les funcions han ser aptes per vectors
    %   a, b: domini de la funcio.
    % OUTPUT:
    %   f: vector solucio.
    %   x: nodes on sha avaluat. Atencio que del primer i del ultim no sen torna el
    %      valor en f.
    [M1, M2, M3, Lhat, x] = crearMatriusODE(n, C, p, q, r);
    mCoef = [C(2, 2), 0;
        0, C(1, 2)];
    esquerra = Lhat + M3 * inv(M2) * mCoef * M1;
    dreta = g(x(2:end - 1)) - M3 * inv(M2) * [C(2, 3); C(1, 3)];

    f = esquerra \ dreta;
