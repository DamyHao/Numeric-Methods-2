%li donem x i y i treu dx/dx, dy/dt
clear all; close all;
exb = @(v)([v(2); -v(1) + 50 * v(2) * (1 - v(1)^2)]);
v0 = [0; 10^ - 2; 10^ - 2];
h = 0.001; % Aquesta es la h aproximada
desiredPoints = 100000;

V1 = ExplicitEulerwTime(v0, h, exb, desiredPoints);

plot(V1(1, :), V1(2, :));
figure;
plot(V1(2, :), V1(3, :));
grid on;
hold on;
figure
V = RK4wTime(v0, h, exb, desiredPoints);
plot(V(1,:), V(2, :), '-k');
hold off;

%% Seccio 2:
close all
%Exacta:
figure;
ff = @(t)((exp(1) / (exp(2) - 1)) .* (exp(t) - exp(-t)));
time = 0:0.001:1;
plot(time, ff(time));
hold on

p = @(x)(0 .* x + 1); % En cas que sigui una funcio constant, cal que apareixi la x en el function handler
q = @(x)(0 .* x);
r = @(x)(0 .* x - 1);
g = @(x)(0.*x);
n = 32;
C = [1, 0, 0; 1, 0, 1];
a = 0;
b = 1;
[f, x] = resoldreODE(n, C, p, q, r, g, a ,b);
plot(x(2: end-1), f, '-ok');
hold off;

%Error: 
semilogy(x(2:end-1), abs(f-ff(x(2:end-1))), '-ok')

%test
ftest = @(x)(2.*x.^2);
[D,x] = chebdiff(64,-4,3);

%{
 figure;
plot(x, ftest(x));
figure
plot(x, D*ftest(x));
figure
plot(x, (D^2)*ftest(x));
%}


