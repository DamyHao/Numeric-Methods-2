function [F, x, lamb] = resoldreODEeig(n, C, p, q,r, a, b)
% Les funcions entregades han de ser les corresponents al domini -1, 1
% Totes les funcions han ser aptes per vectors
[M1,M2,M3,Lhat, x] = crearMatriusODE(n, C, p, q, r, a, b);
mCoef = [ C(2,2), 0 ;
    0, C(1,2)];
esquerra = Lhat + M3*inv(M2)*mCoef*M1;

[EVEC, EVAL] = eig(-esquerra);
lamb = diag(EVAL); % Lamb es un vector perq diag funciona en les dos direccions de conversio.
[foo,ii] = sort(lamb) ; lamb = lamb(ii) ; EVEC = EVEC(:,ii);
eval=diag(lamb);

F = eval*EVEC';

