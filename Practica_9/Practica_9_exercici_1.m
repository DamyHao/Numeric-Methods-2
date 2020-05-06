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
title('Solution found using Newtons method with initial guess u=1')
xlabel('z')
ylabel('u(z)')


%% Section B)
clc 
clear all
n=26;
f0=ones(n-1,1);
U=[];

Norm=[];


n=26;
lambda = exp(1);
p = @(x)(x.*0 +1); 
q =@(x)(0.*x + 0);
r = @(x)(0.*x + 0);
C=[1 0 0; 1 0 0];
a=0; b=1;
normspecial=@(u)(sqrt(cuadratura_cc(a, b, n-2, u.^2)));
mCoef = [ C(2,2), 0 ;  0, C(1,2)];
[M1,M2,M3,Lhat, x] = crearMatriusODE(n, C, p, q, r,a,b);
F = @(f)((Lhat + M3*inv(M2)*mCoef*M1)*f + M3*inv(M2)*[C(2,3);C(1,3)] + lambda.*exp(f));
lambda=[-4:0.1:4];
for i=lambda
    [M1,M2,M3,Lhat, x] = crearMatriusODE(n, C, p, q, r,a,b);
    mCoef = [ C(2,2), 0 ;  0, C(1,2)];
    F= @(f)((Lhat + M3*inv(M2)*mCoef*M1)*f + M3*inv(M2)*[C(2,3);C(1,3)] + i.*exp(f));
    [XK, resd, it] = newtonn(f0, 1e-12, 100, F);
    U=[U XK(:,end)];
    f0=XK(:,end);
    Norm=[Norm norm(XK(:,end))];
end
figure;
plot(lambda,Norm)

%% 

%Manera continuation step

clc;
clear all;
close all;


%Variables to find the firsts solutions
n=26;
eps = 0.1;
lambdas = [-4, -4+eps];
p = @(x)(x.*0 +1); 
q =@(x)(0.*x + 0);
r = @(x)(0.*x + 0);
C=[1 0 0; 1 0 0];
a=0; b=1;
mCoef = [ C(2,2), 0 ;  0, C(1,2)];
[M1,M2,M3,Lhat, x] = crearMatriusODE(n, C, p, q, r,a,b);



aleatoryTimes = 1:1:2;
figure;
S=[];

for lamb = lambdas
    F = @(f)((Lhat + M3*inv(M2)*mCoef*M1)*f + M3*inv(M2)*[C(2,3);C(1,3)] + lamb.*exp(f));
    f0 = ones(n-1, 1);
    [XK, resd, it] = newtonn(f0, 1e-12, 100, F);
    S=[S XK(:,end)];
    plot(x(2:end-1), XK(:, end), '-*');
    hold on
end

title('Exploration at lambda = e using newton with aleatory initial points');
xlabel('z')
ylabel('u(z)')
hold off



%To do the continuation step
s = 1;
itmax = 100;
tol = 1e-10;
eps = 0.1;
y0=[S(:,1); -4];
y1=[S(:,2); -4+eps];
a=0; b=1;
normspecial=@(u)(sqrt(cuadratura_cc(a, b, n-2, u.^2)));
mCoef = [ C(2,2), 0 ;  0, C(1,2)];
[M1,M2,M3,Lhat, x] = crearMatriusODE(n, C, p, q, r,a,b);
funLamb = @(f)((Lhat + M3*inv(M2)*mCoef*M1)*f(1:end-1) + M3*inv(M2)*[C(2,3);C(1,3)] + f(end).*exp(f(1:end-1)));
y=y1;
Y = [];
it = 0;
normspecial=@(u)(sqrt(cuadratura_cc(a, b, n-2, u.^2)));
normes=[];
figure;
while -4<y(end) && y(end)<4 && it < 100
[y, iconv] = continuationStep(funLamb, y0, y1, s, tol, itmax);
y0 = y1;
y1 = y;
Y = [Y, y];
normes=[normes, normspecial(y(1:end-1))];
plot(y(end), norm(y(1:end-1)), 'o');
hold on
it = it +1;
end
hold off

figure;
plot(Y(end, :),normes, 'o');
title('Exploration from lambda=-4 to lambda=4 using secant continuation step');
xlabel('lambda')
ylabel('U')
hold off












