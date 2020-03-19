clear all;
close all;
format long g;

c = (2*9.81)^(-1/2);

funcio = @(x)(((1+(x.^2)./(1-x.^2))./((1-x.^2).^(1/2))).^(1/2));
fun2 = @(x)((1-x.^2).^(-3/4));
interval = sqrt(2)/2;
a = -interval;
b = interval;
Ns = [4,6,8,12];

for n = Ns
       i = cuadratura_cc(a, b, n, funcio);
       i = c*i;
       i2 = cuadratura_cc(0, b, n, fun2);
       i2 = i2*(2/9.81)^(1/2);
end

% b)
Ns = [12,16,20,24];
for n = Ns
       ifej = fejer(-1, 1, n, funcio);
       ifej = c*ifej;
       %i2 = fejer(0, b, n, fun2);
       %i2 = i2*(2/9.81)^(1/2);
end

% c)

const = (2*9.81)^(-1/2);
Cs = [2,4,6];
resultats = [];
for c = Cs
    itan = cuadratura_tanh(-1, 1, c, 16, funcio);
    itan = const*itan;
    resultats = [resultats itan];
end
resultats

% d) Pista : 1-x = z^4
% Funcio sense patologies (no te asímptotes)
funcioSana = @(z)((2-z.^4).^(-3/4));
% Recordar a treure el signe negatiu de dintre l'arrel!!!!
const = (2/9.81)^(1/2);
Cs = [8,12,16,24];
resultats2 = [];
for n = Ns
    itan = cuadratura_cc(0, 1, n, funcioSana);
    itan = 4*const*itan;
    resultats2 = [resultats2 itan];
end

abs(resultats2)

A = [1 1 1; 0 1/2 1; 0 1/4 1];
B = [-1; -1/4; -1/9];
X = linsolve(A, B);

% Problema 2
funci = @(x)(400.*exp(-0.0014.*x).*(cos(1.6*10^(-6).*x.^2)).^2);
xs = [0:1:2000];
plot(xs, funci(xs));

xa = newton2(800, 10^-6, 1000, funci);
a = 0.0014;
H = 400;
b = 1.6*10^-6;
derivadeta = @(x)(-H*exp(-a.*x).*cos(b.*x.^2).*(a.*cos(b.*x.^2)+4*b.*x.*(sin(b*x.^2))));
derivadeta2 = @(x)(a.*cos(b.*x.^2)+4*b.*x.*(sin(b*x.^2)));

xb = newton2(1000, 10^-6, 1000, derivadeta);

novader = @(x)((-H*exp(-a.*x).*cos(b.*x.^2).*(a.*cos(b.*x.^2)+4*b.*x.*(sin(b*x.^2)))).^2);
aintegrar = @(x)(((1+novader(x))./(H-funci(x)).^2));
aintegrar2 = @(x)(((1+((-H*exp(-a.*x).*cos(b.*x.^2).*(a.*cos(b.*x.^2)+4*b.*x.*(sin(b*x.^2)))).^2))./(400-(400.*exp(-0.0014.*x).*(cos(1.6*10^(-6).*x.^2)).^2))).^0.5);

integral = fejer(0, xa, 1000, aintegrar2);
integral2 = anTanh(0, xa, 250, aintegrar2, 1);
%integral3 = cuadratura_tanh(0, xa, 1, 250, aintegrar2)
tt = (2*9.81)^(-1/2)*integral;

funci2 = @(x)((400.*exp(-0.0014.*x).*(cos(1.6*10^(-6).*x.^2)).^2)-20);
I = cuadratura_cc(849.1328, 1143.43, 30, funci2);