function [T,V]=ExplicitEulerTime(fun,t0,h,v0,n)
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

V = zeros(length(v0), n-1);
V(:,1) = v0;
vn = v0;
t = t0;
T=[t];

for i = 1:n-2 % Si li demanem un punt fara 0 iteracions
    a = h * fun(vn);
    vn1 = vn + a;
    V(:, i+1) = vn1; %cada columna es un "punt"
    vn = vn1;
    t = t + h;
    T=[T t];
end

T=t0+[1:n-1].*h;
end