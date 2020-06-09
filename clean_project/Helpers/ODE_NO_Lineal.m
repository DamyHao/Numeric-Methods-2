function F = ODE_NO_Lineal(f)
% Aquest codi et dona les intruccions sobre com crear la funcio que se li 
% ha de ficar al newton per tal de resoldre l'odre no lineal

%Discretitzem la part lineal de l'ode per obtenir terme lineal:
[M1,M2,M3,Lhat, x] = crearMatriusODE(n, C, p, q, r,a,b);
mCoef = [ C(2,2), 0 ;
    0, C(1,2)];
%terme lineal:
lineal = (Lhat + M3*inv(M2)*mCoef*M1)*f +M3*inv(M2)*[C(2,3);C(1,3)];

%TERME NO LINEAL (MODIFICAR EN CADA CAS)
nolineal= f.*D*f;

%Funcio pel Newton
F=lineal-nolineal;  %VIGILAR AMB ELS SIGNES


%% EXEMPLE D'ÚS (exam final 2019 exercici 2):

lamb=18;
C=[1,0,0;...
   1,0,0];
n=26; a=0; b=1;
p=@(x)(x.*0+1/lamb);
q=@(x)(x.*0);
r=@(x)(x.*0 + 1);

[M1,M2,M3,Lhat, x] = crearMatriusODE(n, C, p, q, r,a,b);
mCoef = [ C(2,2), 0 ;
    0, C(1,2)];
%lineal = (Lhat + M3*inv(M2)*mCoef*M1)*f +M3*inv(M2)*[C(2,3);C(1,3)];
[D,x] = chebdiff(n,a,b);
%nolineal=f.*D*f;
%newtonfunction=@(f)(lineal-nolineal);
newtonfunction=@(f)((Lhat + (M3*inv(M2))*mCoef*M1)*f +M3*inv(M2)*[C(2,3);C(1,3)]-f.*D(2:n,2:n)*f); %important D(2:n,2:n)


x0=sin(pi.*x(2:n));

[XK, resd, it] = newtonn(x0, 1-8, 100, newtonfunction);
plot(x(2:n),XK(:,end))


