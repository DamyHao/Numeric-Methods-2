%% Sense Freganent
clear all;
close all;
Vzo = 0;
zo = 120000;

% Definim les funcions:
znext = @(Zn, Vn, t, a)(Zn + Vn*t + 0.5*a*(t^2));
Vznext = @(Vn, t, a)(Vn + a*t);

it = 0;
a = -9.81;
inTemps = 0.01;
v0 = Vznext(Vzo, 0, a);
z0 = znext(zo, v0, 0, a);
t = 0;

velocitats = [];%zeros(15641);
zetas = [];%zeros(15641);

 while z0 > 0 && it < 1000000
     t = t + inTemps;
     it = it + 1;
     v0 = Vznext(v0, inTemps, a);
   	 z0 = znext(z0, v0, inTemps, a);
     zetas = [zetas z0];
     velocitats = [velocitats v0];
 end 
 
 temps = [0:inTemps:t-inTemps];
 plot(temps, zetas);
 hold on;
 plot(temps, velocitats);
 hold off;
 
 % 3
 % Utilitzant el temps del apartat anterior.
 
 Ynext = @(y, t, v)(y + v*t);
 ys = [];
 yActual = 0;
 Vyo = 7833.3;
 for t = temps
     yActual = Ynext(yActual, inTemps, Vyo);
     ys = [ys yActual];
 end
%  while z0 > 0 && it < 1000000
%      t = t + inTemps;
%      it = it + 1;
%      v0 = Vznext(v0, inTemps, a);
%    	 z0 = znext(z0, v0, inTemps, a);
%      zetas = [zetas z0];
%      velocitats = [velocitats v0];
%  end 
figure;
plot(temps, ys);

figure;
plot(ys, zetas);

%% Apartat 2:
clear all;

znext = @(Zn, Vn, t, a)(Zn + Vn*t + 0.5*a*(t^2));
Vznext = @(Vn, t, a)(Vn + a*t);

Voz = 0;
g = -9.81;
an = g;
massa = 69103;
incTemps = 0.001;

xn1 = znext(120000, Voz, incTemps, an);
vn1 = Vznext(Voz, incTemps, an);

b = 115;
F(v) = @(v)(115*v);

an =  (g * massa + F(vn1))/massa;

