function [x1,x2] = arrels(a,b,c);
d=b^2-4*a*c;
if a~=0
    x1 = (-b+sqrt(d))/(2*a);
    x2 = (-b-sqrt(d))/(2*a);
elseif b~=0
    x1 = -c/b;
    x2 = 'no existeix'
    disp('equació lineal amb arrel única')
else c~=0
    disp('equació impossible');
end
