close all;
clc;
clear all;

figure
p = @(x)(0 .* x + 1); % En cas que sigui una funcio constant, cal que apareixi la x en el function handler
q = @(x)(0 .* x + 1);
r = @(x)(0 .* x + 100);
g = @(x)(1 .* x);
n = 1280;
CC = [1, 0, 0; 1, 0, 0];
a = -1; b = 1;
[f, x] = resoldreODE(n, CC, p, q, r, g, a, b);
plot(x(2:end-1), f, '-k.');

figure
p = @(x)(0 .* x + 1); % En cas que sigui una funcio constant, cal que apareixi la x en el function handler
q = @(x)(0 .* x + 1);
r = @(x)(0 .* x + 100);
g = @(x)(1 .* x);
n = 1280;
CC = [0, 1, 0; 0, 1, 0];
a = -1; b = 1;
[f, x] = resoldreODE(n, CC, p, q, r, g, a, b);
plot(x(2:end-1), f, '-k.');

%% Exercici 2:
close all
format long g;

v0 = [-2.5; -2.5];

[XK, resd, it] = newtonn(v0, 1e-6, 150, @newtonFunction);
XK(:, end)

v0 = [-2.5; 2.5];
% C
[XK2, resd, it] = newtonn(v0, 1e-10, 150, @newtonFunction);
XK2(:, end)

v0 = [-2.5; 0];

[XK3, resd, it] = newtonn(v0, 1e-9, 150, @newtonFunction);
XK3(:, end)


h = 0.01;
vn0 = [0;-2.8 ; 0 ; 3.1 ;0];
finalTime = 60;
   %   Tfinal sera vn(0) + h*(desiredPoints-1)
desiredPoints = finalTime/h + 2;
V = RK4wTime(vn0, h, @rkFun, desiredPoints);

figure;
vn0 = [0;-2.8 ; 0 ; 3.1 ;0];
finalTime = 5;
desiredPoints = finalTime/h + 1;
V2 = RK4wTime(vn0, h, @rkFun, desiredPoints);
plot(V2(1,:), V2(2,:));

figure;

xt = V(2,1:end)';

N = length(xt);
%xt = sin(1:N);
fourierSpace = dftmat(xt);
k = -N/2:(N/2-1);
wk = 2*pi/(N*h).*k; %spectrum of frequencies
semilogy(wk, abs(fourierSpace));


XK2(:, end);
xEqui = [XK(1,end), 0 , XK(2, end), 0];
xEqui2 = [XK2(1,end); 0 ; XK2(2, end); 0];
DF = jaco(@rkFun, xEqui);
DF2 = jaco(@rkFun, xEqui2);
[VEC, VAL] = eig(DF2);
disp(abs(diag(VAL)))