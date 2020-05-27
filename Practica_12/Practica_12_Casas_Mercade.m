%% Practica_20_Casas_Mercade

%% Section A

close all 
clear all
clc

N=128;
haches=[0.1, 0.01];
D1=dftdiffmat1(N);
D2=dftdiffmat2(N);
Burger=@(f)(0.1*D2*f+f.*D1*f);
for h=haches
    x=(2*pi/N)*[0:N-1]';
    f0=exp(-10.*(cos(x/2)).^2);
    f = ExplicitEuler(f0, h,Burger ,100);
end

figure;
%The stability region of the Explicit Euler (RK1) is a circle of radius 1
%centered at z=-1
pos = [-2 -1 2 2]; %Com si fos un rectangle
r = rectangle('Position',pos,'Curvature',[1 1]);
r.FaceColor = [0.752, 0.984, 0.674];
text(-1.5,0.5,  ['Stable region'], 'FontSize', 12)
title('Explicit Euler')
axis([-5 0.5 -2 2])
grid on
hold on

%Now we compute the jacobian evaluated at the initial condition f0
DF = jaco(Burger, f0);
[evec,eval]=eig(DF);
z=diag(eval).*haches;

plot(z,'*')
hold on

% This integrations is not stable neither for h=0.1 nor h=0.01 because
% in both cases there are points outside the stability region of euler's
% explicit method.

%% Section B)

h = 0.1;
fBurger=@(t,f)(0.1*D2*f+f.*D1*f);


figure;
%The stability region of the Explicit Euler (RK1) is a complex circle of radius 1
% centered at z=1
pos = [0 -1 2 2]; %Com si fos un rectangle
r = rectangle('Position',pos,'Curvature',[1 1]);
r.FaceColor = [1, 0.2, 0.274];
text(0.4,0,  ['Unstable region'], 'FontSize', 12)
title('Implicit Euler')
grid on
axis([-4 2.5 -2 2])
hold on

%Now we compute the jacobian evaluated at the initial condition f0
DF = jaco(Burger, f0);
[evec,eval]=eig(DF);
z=diag(eval).*haches;
plot(z,'*')
hold on


figure;
times = [2,4];
h = 0.1;
Results=[];
for t = times
    steps = t/h;
    [time, result] = bdf1(fBurger, 0, h, f0, steps);
    Results=[Results, result(:,end)];
end

plot(x, Results);
legend('t=2','t=4')
xlabel('x')
ylabel('y')





