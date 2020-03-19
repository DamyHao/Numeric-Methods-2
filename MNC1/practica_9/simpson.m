function I = simpson(a, b, m, fun)
%   Retorna la integral definida utilitzant el metode de integracio compost.
%   a i b son els intervals d'integracio.
%   fun es la funcio que volem integrar.

H=(b-a)/m ; %distancia entre nodes.
k=[0:2*m];
xk=[a+(H.*k)./2];

% Sumem 1 al primer terme de cada sumatori perque les arrays comencen a 1.
sum1 = 0;
sum2 = 0;
s = 0;
for r=3:2:2*m
    sum1= sum1 + fun(xk(r));
    s = r-1;
    sum2 = sum2 + fun(xk(s));
end
sum2 = sum2 + fun(xk(s+2));
%EXPLICACIO

I = (H/6) * (fun(xk(1)) + 2*sum1 + 4*sum2 + fun(xk(end)));


