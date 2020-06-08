function v = newtonFunction (in)
v = zeros(2, 1);
x = in(1);
y = in(2);

v(1) =  4*(x^2 + y - 11)*x + 2*(x + y^2 -7);
v(2) =  2*(x^2 + y - 11) + 4*(x + y^2 -7)*y;

end