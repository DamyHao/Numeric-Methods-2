%% 4) Proves
clear all;
x=0; y=yo-h; z=zo;
dV = @(th)(-(dens*radi)./(sqrt((radi^2).*(cos(th)).^2 + (y - radi.*sin(th)).^2 + z.^2)));
Vy2 = -(cuadratura_cc(0, 2*pi, 20, dV);

Fy = -(Vy1-Vy2)/(2*h);

x=0; y=yo; z=zo+h;
dV = @(th)(-(dens*radi)./(sqrt((radi^2).*(cos(th)).^2 + (y - radi.*sin(th)).^2 + z.^2)));
Vz1 = -(cuadratura_cc(0, 2*pi, 20, dV);

x=0; y=yo; z=zo-h;
dV = @(th)(-(dens*radi)./(sqrt((radi^2).*(cos(th)).^2 + (y - radi.*sin(th)).^2 + z.^2)));
Vz2 = -(cuadratura_cc(0, 2*pi, 20, dV);

Fz = -(Vz1-Vz2)/(2*h);

epsi = 0.001;
nextPoint = @(yo,z)([yo+epsi*])

