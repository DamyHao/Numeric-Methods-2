%% Exercici_1_BVP_Casas_Mercade
% Non linear differential operator
%% Section A)
clc
clear all
close all

n=26;
f0=ones(n-1,1);
[XK, resd, it] = newtonn(f0, 1e-12, 100, @newtonFunction);

[D, x]= chebdiff(n, 0 ,1);
plot(x(2:end-1),XK(:,end));
a=0; b=1;
normSpecial=@(u)(sqrt(cuadratura_cc(a, b, n-2, u.^2)));

disp(normSpecial(XK(:,end)))
title('Solution found using Newtons method with initial guess u=1')
xlabel('z')
ylabel('u(z)')


%% Section B)
clc
clear all


%per crear la funcio
n=26; a=0; b=1;
p = @(x)(x.*0 +1);
q =@(x)(0.*x + 0);
r = @(x)(0.*x + 0);
C=[1 0 0; 1 0 0];
mCoef = [ C(2,2), 0 ;  0, C(1,2)];
[M1,M2,M3,Lhat, x] = crearMatriusODE(n, C, p, q, r,a,b);

%norma U
normSpecial=@(u)(sqrt(cuadratura_cc(a, b, n-2, u.^2)));

U=[];
Norm=[];
f0=ones(n-1,1);
initialLamb = -4;
eps = 0.1;
lambdas=[initialLamb, (initialLamb + eps)];
initialGuesss = [];

% Trobem les dos condicions inicials
for la=lambdas
    F = @(f)((Lhat + M3*inv(M2)*mCoef*M1)*f + M3*inv(M2)*[C(2,3);C(1,3)] + la.*exp(f));
    [XK, resd, it] = newtonn(f0, 1e-10, 100, F);
    initialGuesss = [initialGuesss, [XK(:, end); la]];
end


%To do the continuation step
s = 1; itmax = 500; tol = 1e-10;
y0=initialGuesss(:,1);
y1=initialGuesss(:,2);

%normes=[normSpecial(y0(1:end-1)), normSpecial(y1(1:end-1))];
%lamb=[y0(end), y1(end)];
%Y = [y0, y1];


funLamb = @(f)((Lhat + M3*inv(M2)*mCoef*M1)*f(1:end-1) + M3*inv(M2)*[C(2,3);C(1,3)] + f(end).*exp(f(1:end-1)));
iterator = 0;
y= y0;
Y = [];

figure;
while -5<y(end) && y(end)<4 && iterator < itmax
    [y, iconv] = continuationStep(funLamb, y0, y1, s, tol, itmax);
    y0 = y1;
    y1 = y;
    Y = [Y, y];
    iterator = iterator +1;
    yy =  normSpecial(y(1:end-1));
    plot(y(end), yy, 'o');
    hold on;
end
hold off

% plot(lamb,normes, 'o');
% title('Exploration from lambda=-4 to lambda=4 using secant continuation step');
% xlabel('lambda')
% ylabel('U')
% hold off












