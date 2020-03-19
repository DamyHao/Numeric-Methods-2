function [y1,z1] = nextPoint(yo,zo,xk)
%NEXTPOINT Calcula el conjunto de puntos equipotenciales.
    % A partir de un punto concreto (yo,zo) cualquiera sometido
    % un potencial determinado, calculamos el siguiente punto que se
    % haya bajo el mismo potencial.

%Datos de los que disponemos: 
th = xk; x = 0;
dens = 1; radi = 1;
a=0; b=2*pi;
N=20;
dV = @(th,dens,radi,x,y,z)(-(dens*radi)./(sqrt((radi^2).*(cos(th)).^2 + (y - radi.*sin(th)).^2 + z^2)));

%Cálculos para obtener el valor de las diferencias finitas:
h=0.01;

th = xk; y=yo+h; z=zo;
fx = dV(th,dens,radi,x,y,z);
V1 = clenshawcurtis_p10(a, b, N, fx);
th = xk; y=yo-h; z=zo;
fx = dV(th,dens,radi,x,y,z);
V2 = clenshawcurtis_p10(a, b, N, fx);

Fy = -(V1-V2)/2*h;

th = xk; y=yo; z=zo+h;
fx = dV(th,dens,radi,x,y,z);
V1 = clenshawcurtis_p10(a, b, N, fx);
th = xk; y=yo; z=zo-h;
fx = dV(th,dens,radi,x,y,z);
V2 = clenshawcurtis_p10(a, b, N, fx);

Fz = -(V1-V2)/2*h;

epsilon = 0.001;
f_uni_y = (-Fz)/sqrt(Fy^2 + Fz^2);
f_uni_z = Fy/sqrt(Fy^2 + Fz^2);

y1 = yo + epsilon.*f_uni_y;
z1 = zo + epsilon.*f_uni_z;
end

