function [M1,M2,M3,Lhat, x] = crearMatriusODE(n, C, p, q, r,a,b)
% Crea les matrius necessaries per les funcions resoldreODE i resoldreODEeig estandard
% (de la carpeta PRACTICA_8_NO). Per veure mes informacio, ves a aquestes funcions.
% S'hi s'han de resoldre ODEs que no tenen totes les boundary conditions,
% anar a PRACTICA8. (No es trobara en clean_project)

% Fem la diferenciacio Chebyshev, que ens retorna també els nodes:
[D,x] = chebdiff(n,a,b);
P = diag(p(x));
Q=diag(q(x));
R=diag(r(x));
L=P*D^2+Q*D+R;

Lhat=L(2:end-1,2:end-1);

M1=-[D(1,2:end-1); D(end,2:end-1)];

% Factor que fa que tot el canvi de a,b esta amagat en D. No cal preucuparsen
M2 = [C(2,1) + C(2,2)*D(1,1), C(2,2)*D(1,end);
      C(1,2)*D(end, 1), C(1,1) + C(1,2)*D(end,end)];
    
M3= [L(2:end-1,1), L(2:end-1,end)];

end