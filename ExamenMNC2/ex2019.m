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
plot(V(1, :), V(2, :), '-k');
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
g = @(x)(0 .* x);
nn = 32;
CC = [1, 0, 0; 1, 0, 1];
a = 0;
b = 1;
[f, xx] = resoldreODE(nn, CC, p, q, r, g, a, b);
plot(xx(2:end - 1), f, '-ok');
hold off;

%Error:
semilogy(xx(2:end - 1), abs(f - ff(xx(2:end - 1))), '-ok')

%test

%{
ftest = @(x)(2 .* x.^2);
[D, x] = chebdiff(64, -4, 3);
%}

%{
figure;
plot(x, ftest(x));
figure
plot(x, D * ftest(x));
figure
plot(x, (D^2) * ftest(x));
%}
%%
close all;
clc;
clear all;
lambdas = [18, 16, 14];
global n D D2 I lamb x C
n = 64;
a = 0;
b = 1;
[D, x] = chebdiff(n, a, b);
%D2 = D^2;

C = [1, 0, 0; 1, 0, 0];

figure

for lamb = lambdas
    f0 = (sin(pi * x(2:end - 1)));
    [XK, resd, it] = newtonn(f0, 1e-10, 100, @fap);
    plot(x(2:end - 1), XK(:, end));
    [foo, ii] = max(XK(:, end));
    text(x(ii), XK(ii, end), ['lambda = ', num2str(lamb)]);
    hold on;
end

hold off

% Ultim apartat Trobem inicials:
lamb = 18;
f0 = (sin(pi * x(2:end - 1)));
[XK, resd, it] = newtonn(f0, 1e-10, 100, @fap);
y0 = [XK(:, end); lamb];

lamb = 17.9;
f0 = (sin(pi * x(2:end - 1)));
[XK, resd, it] = newtonn(f0, 1e-10, 100, @fap);
y1 = [XK(:, end); lamb];

s = 1;
maxIt = 1000;
iterations = 0;
Y = [];

while s > 0 && iterations <= maxIt
    [y, iconv] = continuationStep(@fapLambda, y0, y1, s, 1e-6, 100);

    if iconv == 1% No hem aconseguit solució i ajustem s
        s = s - 0.1; % Si la s arriba a 0 desistirem i no buscarem mes solucions
    else
        y0 = y1;
        y1 = y;
        Y = [Y, y]; %solucions
    end
    plot(y(end), max(y(1:end-1)), 'o');
    hold on;
    iterations = iterations + 1;
end


%{
 plot(Y(end, :), Y(1:2, :), 'o');
title('Continuation step');
xlabel('Alpha')
ylabel('Phi')
hold on
 
%}

% Per solucionar 1/lamb * f'' + f + f*f' = 0 separem en la part lineal i no lineal
% Per la lineal farem servir operador
function F = fap(f)
    global n D I lamb x C
    % Fem la diferenciacio Chebyshev, que ens retorna també els nodes:
    p = @(xx)(0 .* xx + 1 / lamb); % En cas que sigui una funcio constant, cal que apareixi la x en el function handler
    q = @(xx)(0 .* xx);
    r = @(xx)(0 .* xx + 1);
    P = diag(p(x));
    Q = diag(q(x));
    R = diag(r(x));
    L = P * D^2 + Q * D + R; % Operador
    % Les constants del canvi de variable es posen automaticament
    % No cal fer L = 4*P * D^2 + 2 *Q * D + R; % Operador
    Lhat = L(2:end - 1, 2:end - 1);

    M1 = -[D(1, 2:end - 1); D(end, 2:end - 1)];

    % Factor que fa que tot el canvi de a,b esta amagat en D. No cal preucuparsen
    M2 = [C(2, 1) + C(2, 2) * D(1, 1), C(2, 2) * D(1, end);
        C(1, 2) * D(end, 1), C(1, 1) + C(1, 2) * D(end, end)];

    M3 = [L(2:end - 1, 1), L(2:end - 1, end)];
    mCoef = [C(2, 2), 0; 0, C(1, 2)];
    esquerra = Lhat + M3 * inv(M2) * mCoef * M1; % Si gamma no fos 0 en el anunciat opino que sortiria un altre terme (el de la M3*M2^-1)
    % nonLineal = 2 * f .* (D(2:end-1, 2:end-1) * f);
    nonLineal = f .* (D(2:end - 1, 2:end - 1) * f);
    F = esquerra * f - nonLineal;
end

% Funcio amb lambda per el ultim apartat
function F = fapLambda(f)
    global n D I x C
    % Fem la diferenciacio Chebyshev, que ens retorna també els nodes:
    p = @(xx)(0 .* xx + 1 / f(end)); % En cas que sigui una funcio constant, cal que apareixi la x en el function handler
    q = @(xx)(0 .* xx);
    r = @(xx)(0 .* xx + 1);
    P = diag(p(x));
    Q = diag(q(x));
    R = diag(r(x));
    L = P * D^2 + Q * D + R; % Operador
    % Les constants del canvi de variable es posen automaticament
    % No cal fer L = 4*P * D^2 + 2 *Q * D + R; % Operador
    Lhat = L(2:end - 1, 2:end - 1);

    M1 = -[D(1, 2:end - 1); D(end, 2:end - 1)];

    % Factor que fa que tot el canvi de a,b esta amagat en D. No cal preucuparsen
    M2 = [C(2, 1) + C(2, 2) * D(1, 1), C(2, 2) * D(1, end);
        C(1, 2) * D(end, 1), C(1, 1) + C(1, 2) * D(end, end)];

    M3 = [L(2:end - 1, 1), L(2:end - 1, end)];
    mCoef = [C(2, 2), 0; 0, C(1, 2)];
    esquerra = Lhat + M3 * inv(M2) * mCoef * M1; % Si gamma no fos 0 en el anunciat opino que sortiria un altre terme (el de la M3*M2^-1)
    % nonLineal = 2 * f .* (D(2:end-1, 2:end-1) * f);
    nonLineal = f(1:end-1) .* (D(2:end - 1, 2:end - 1) * f(1:end-1));
    F = esquerra * f(1:end-1) - nonLineal;
end
