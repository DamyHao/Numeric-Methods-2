function y = Afun (x)

%funció que et permet trobar Ax sense haver de calcular la matriu A,

%{
fiquem un exemple:
%}

R=3;V=1; n=30; Ax=[]; Ax = [Ax ; 2*R*x(1) + R*x(2)]; Ax = [Ax ; x(1)-x(2)-x(3)] ;

for j = 3:n
    Ax = [Ax ; 2*R*x(2*j-3) + R*x(2*j-2) - R*x(2*j-4)] ;
    Ax = [Ax ; x(2*j-3) - x(2*j-2) - x(2*j-1)] ;
end



end