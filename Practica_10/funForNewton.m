function r = funForNewton(theta)
%Funcio per tirar newton en una variable a la practica 10:
disp(theta);
h = 0.0001;
points = 2/h;
initial = [ 0; 0; sqrt(2).*cos(theta); sqrt(2).*sin(theta)];
sol = RK4(initial, h, @gravFunctionV, points);
rs = sol(1:2, end);
r = norm(rs);
end