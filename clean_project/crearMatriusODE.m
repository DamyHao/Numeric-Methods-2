function [M1,M2,M3,Lhat, x] = crearMatriusODE(n, C, p, q, r,a,b)
% Input: n: nombre de nodes Chebyshev,
% C: Matriu 2*3 amb les Robin conditions
% Funcions p, q, r:  En cas que sigui una funcio constant, cal que apareixi la x en el function handler

% Fem la diferenciacio Chebyshev, que ens retorna també els nodes:
[D,x] = chebdiff(n,a,b);
P = diag(p(x));
Q=diag(q(x));
R=diag(r(x));
L=P*D^2+Q*D+R;

Lhat=L(2:end-1,2:end-1);

M1=-[D(1,2:end-1);D(end,2:end-1)];

M2 = [C(2,1) + C(2,2)*D(1,1), C(2,2)*D(1,end);
      C(1,2)*D(end, 1), C(1,1) + C(1,2)*D(end,end)];
    
M3= [L(2:end-1,1), L(2:end-1,end)];


end