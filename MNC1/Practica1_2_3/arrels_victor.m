function [x1,x2] = arrels_victor(a, b, c)
discriminant = b^2-4*a*c;

if a = 0
    disp('equació lineal; única solució');
    x1 = (-c)/b;
elseif b = 0
    x1 = sqrt((-c)/a);
    x2 = -sqrt(-c/a);
elseif c = 0
    disp('equació lineal; única solució');
    x1 = (-c)/a;
elseif a=0; b=0;
    disp('equació impossible');
elseif a = 0; b = 0;
    x1 = 0;
else
    x1 = (-b+sqrt(discriminant))/(2*a);
    x2 = (-b-sqrt(discriminant))/(2*a);
end