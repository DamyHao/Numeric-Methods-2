function [F, x, lamb] = resoldreODEeig(n, C, p, q,r, a, b)
% Resol una equacio diferencial de tipus boundary en el interval a,b.
% ATENCIO: fiquem les C, p, q, r, al domini de (a,b)
% INPUT:
%   n: nombre de nodes
%   C: matriu dels coeficients C
%       {C_11*f(a) + C_12*f(a) = C_13*f(a)
%       {C_21*f(b) + C_22*f(b) = C_23*f(b)
%   p, q, r:funcions del boundary problem p(x)f'' + q(x)f' + r(x)f = lambda*f
%           En cas que sigui una funcio constant, cal que apareixi la x en 
%           el function handler. Les funcions entregades han de ser les 
%           corresponents al domini a, b.
%           Totes les funcions han ser aptes per vectors
%   a, b: domini de la funcio.
% OUTPUT:
%   F: Matriu de eigenfunctions. Cada eigenfunction es una fila.
%      Estaran ordenades per ordre ascendent del seu eigenvalue associat.
%       Exemple: plot(x(2:end-1), 4*F(1:3,:)); (3 primeres)
%   x: nodes on sha avaluat. Atencio que del primer i del ultim no sen torna el
%      valor en F.
%   lamb: vector de eigenvalues ordenats ascendentment

[M1,M2,M3,Lhat, x] = crearMatriusODE(n, C, p, q, r, a, b);
mCoef = [ C(2,2), 0 ;
    0, C(1,2)];
esquerra = Lhat + M3*inv(M2)*mCoef*M1;

[EVEC, EVAL] = eig(-esquerra);
lamb = diag(EVAL); % Lamb es un vector perq diag funciona en les dos direccions de conversio.
[foo,ii] = sort(lamb) ; lamb = lamb(ii) ; EVEC = EVEC(:,ii);

F = EVEC';

