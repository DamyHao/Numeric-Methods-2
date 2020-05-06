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
normspecial=@(u)(sqrt(cuadratura_cc(a, b, n-2, u.^2)));

disp(normspecial(XK(:,end)))
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
normspecial=@(u)(sqrt(cuadratura_cc(a, b, n-2, u.^2)));

U=[];
Norm=[];
f0=ones(n-1,1);
lambda=[];

for i=[-4:0.1:4]
    F = @(f)((Lhat + M3*inv(M2)*mCoef*M1)*f + M3*inv(M2)*[C(2,3);C(1,3)] + i.*exp(f));
    [XK, resd, it] = newtonn(f0, 1e-10, 100, F);
    if XK(:,end)<5
    U=[U XK(:,end)];
    f0=XK(:,end);
    Norm=[Norm normspecial(XK(:,end))];
    lambda=[lambda i];
    end
end
figure;
plot(lambda,Norm,'o')




%To do the continuation step
s = 1; itmax = 500; tol = 1e-13;
y0=[U(:,1); lambda(1)];
y1=[U(:,2); lambda(2)];

normes=[normspecial(y0(1:end-1)), normspecial(y1(1:end-1))];
lamb=[y0(end), y1(end)];
Y = [y0, y1];

funLamb = @(f)((Lhat + M3*inv(M2)*mCoef*M1)*f(1:end-1) + M3*inv(M2)*[C(2,3);C(1,3)] + f(end).*exp(f(1:end-1)));

y=y1;
it = 0;


figure;
while -4<y(end) && y(end)<4 && it < itmax
    
[y, iconv] = continuationStep(funLamb, y0, y1, s, tol, itmax);
y0 = y1;
y1 = y;
Y = [Y, y];

normes=[normes, normspecial(Y(1:end-1,end))];
lamb=[lamb y(end)];

hold on
it = it +1;
end
hold off

figure;
plot(lamb,normes, 'o');
title('Exploration from lambda=-4 to lambda=4 using secant continuation step');
xlabel('lambda')
ylabel('U')
hold off












