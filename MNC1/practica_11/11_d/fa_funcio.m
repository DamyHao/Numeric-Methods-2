function error = fa_funcio(theta)
% Aproximación del periodo del péndulo para pequeñas oscilaciones:
T_posc = 2*pi*((1/9.81)^(1/2));

k = sin(theta/2);
K = @(phi)(1./((1-(k.^2).*(sin(phi)).^2)).^(1/2));
a = 0;
b = pi/2;
% m son els nodes, posarem uns 200 per probar.
m = 100;
T_real = 4*((1/9.81)^(1/2))*fejer(a, b, m, K);


error = ((abs(T_real - T_posc)/T_real) - 0.05);
end

