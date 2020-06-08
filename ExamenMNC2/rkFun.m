function out = rkFun(in)
out = zeros(4,1);
x = in(1);
y = in(3);
px = in(2);
py = in(4);

%dx / dt:
out(1) = px;
% dpx/dt
out(2) = -(4*(x^2 + y - 11).*x + 2*(x + y^2 -7));
out(3) = py;
out(4) = -(2*(x.^2 + y - 11) + 4*(x + y^2 -7).*y);

end
