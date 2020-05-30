function [T,X]=ExplicitEulerTime(fun,t0,h,x0,n)
%INPUTS:
%   fun=funció a integrar 
%   t0= temps d'inici de la integració
%   h=time-step
%   x0=punt inicial columna
%   n=number of points (n-1 steps) inclos el putn inciial x0
%   so the output is (x0;x1;x2...;xn)

%OUTPUTS:
%   T=vector de temps
%   X=matriu on cada columna és un punt

%COMMENT: si volem obtenir el mateix punt final, és a dir el mateix T
%final, per haches diferents, tambe shan de modificar el nombre de punts,
%per exemple si per h=0.1 tenim n=500, si ara canviem a h=0.01 hem de
%canviar a n=5000

X = zeros(length(x0), n-1);
X(:,1) = x0;
xn = x0;
TT=[1];

for i = 1:n-2 % Si li demanem un punt fara 0 iteracions
    a = h * fun(xn);
    xn1 = xn + a;
    X(:, i+1) = xn1; %cada columna es un "punt"
    xn = xn1;
end

T=t0+[1:n-1].*h;
end