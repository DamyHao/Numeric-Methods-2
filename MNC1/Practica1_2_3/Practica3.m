%% P3 Casas Victor
%% 11
function superficie = superficie_cilindre (radi, altura);
area_base = pi*(radi^2);
circumferencia = 2*pi*radi;
superficie = 2*(area_base)*circumferencia*altura;
end

function volum = volum_cilindre(radi, altura);
volum = pi*(radi^2)*altura;
end

%% 12
function [x1,x2] = arrels_victor(a, b, c)
discriminant = b^2-4*a*c;

if a == 0
    disp('equació lineal; única solució');
    x1 = (-c)/b;
elseif b == 0
    x1 = sqrt((-c)/a);
    x2 = -sqrt(-c/a);
elseif c == 0
    disp('equació lineal; única solució');
    x1 = (-c)/a;
elseif a==0; b==0;
    disp('equació impossible');
elseif a == 0; b == 0;
    x1 == 0;
else
    x1 = (-b+sqrt(discriminant))/(2*a);
    x2 = (-b-sqrt(discriminant))/(2*a);
end
end

%%13
%% a)
function out = fibonacci(n)
if n <= 0
    out = 0;
elseif n == 1 || n == 2
    out = 1;
else
    out = fibonacci(n-1) + fibonacci(n-2);
end
end

%% b)
function [num_Au, iter] = numAu(tol, max_iterations)
n1 = 1;
n2 = 1;
it = 1;
result1 = 1;
result2 = 1;
aux = 0;
continuar = 0;

while it < max_iterations && continuar == 0
    aux = n1 + n2;
    n1 = n2;
    n2 = aux;
    result2 = n2/n1;
    
    %: %disp('Resultat 1:')
    %: %disp(result1);
    %: %disp('Resultat 2:')
    %: %disp(result2);
    %: %Comprovar si estem a tolerancia m�nima. 
    if abs(result2 - result1) < tol
        continuar = 1;
    end
         
    %: %Actualitzar comptadors:
    it = it + 1;
    result1 = result2;
end
iter = it;
num_Au = result2;
end

%% 14
function out = derivadaEnPunt(F, x)
%: %Amb els m�s petits que menys 16 no funciona (surt 0).
%: %Amb alguns nombres entre -10 y -16 te un comportament extrany ja que
%: %sembla tenir menys precisi�.
minNumber = 10^-10;
out = (F(x + minNumber) - F(x)) / minNumber;

end