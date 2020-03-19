%% P10: CASAS Y JIMENEZ
%% 1) Potencial en el eje:
clear all;
close all;
format long g;
N=20;
%dV = @(th, y, z, alfa, radi, dens)(dens*radi)/(sqrt((radi^2)*(cos(th))^2 
%+ (y - radi*sin(th))^2 + z^2));

radi=1;
dens=1;
z=0;
y=0;
dV = @(th)(-(dens*radi)./(sqrt((radi^2).*(cos(th)).^2 + (y - radi.*sin(th)).^2 + z^2)));

V = cuadratura_cc(0, 2*pi, 20, dV);

%% 2) Potencial en el plano:
x = 0;
Ys = [0.25 0.75 1.5];
figure;
for y = Ys
    V=[];
    Zs = [-5:0.1:5];
    for z = Zs
        dV = @(th)(-(dens*radi)./(sqrt((radi^2).*(cos(th)).^2 + (y - radi.*sin(th)).^2 + z.^2)));
        V = [V cuadratura_cc(0, 2*pi, 20, dV)];
    end
    plot(Zs,V);
    hold on;
end
hold off;

%% 3) Curvas equipotenciales sobre el plano:
clear all;
clc

xo = 0;

h=0.01;
% Definim la funcio diferencial:
radi=1;
dens=1;
dV = @(th, y, z)(dens*radi)./(sqrt((radi^2).*(cos(th)).^2 + (y - radi.*sin(th)).^2 + z^2));

% Creem els nodes de Chebyshev (per N=20):
N = 20;
a = 0;
b = 2*pi;
j = [0:1:N];
xcheb = cos(j.*pi./N);
ths = a + ((b-a)./2).*(xcheb+1); % Perque aquest +1??????

%dV = @(th)(-(dens*radi)./(sqrt((radi^2).*(cos(th)).^2 + (y - radi.*sin(th)).^2 + z.^2)));
% En treiem un vector de Chebyshev:

figure;

for pp = 0.1:0.1:0.3
    yo = pp;
    zo = pp;
    mida = pp*70000;
    % optimitzacions per fer que el programa vagi més rapid:
    Ys = [];
    Zs = [];
    
    for i = 0:1:mida
        xkascc = dV(ths, yo+h, zo);
        Vy1 = -(cccVector(0, 2*pi, xkascc));
        
        xkascc = dV(ths, yo-h, zo);
        Vy2 = -(cccVector(0, 2*pi, xkascc));
        Fy = -(Vy1 -Vy2)/(2*h);
        
        xkascc = dV(ths, yo, zo+h);
        Vz1 = -(cccVector(0, 2*pi, xkascc));
        
        xkascc = dV(ths, yo, zo-h);
        Vz2 = -(cccVector(0, 2*pi, xkascc));
        Fz = -(Vz1 -Vz2)/(2*h);
        
        punt = [yo, zo];
        punt2 = nextPoint(punt, Fy, Fz);
        % Finalment afegim el punt equipotencial trobat als nostres vectors de
        % punts:
        yo = punt2(1);
        zo = punt2(2);
        %posibles optimitzacions:
        %Ys(i+1) = yo;
        %Zs(i+1) = zo;
        Ys = [Ys yo];
        Zs = [Zs zo];
        
    end
    plot(Ys, Zs);
    hold on;
end

hold off;

disp("Acabat!")