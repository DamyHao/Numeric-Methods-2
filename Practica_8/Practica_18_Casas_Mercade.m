% clear all;
% clc;
% close all;
% 
% C = [1, -1, 1 + pi;
%     1, 1, 1-pi];
% 
% p = @(x)(0.*x + 1); % En cas que sigui una funcio constant, cal que apareixi la x en el function handler
% q =@(x)(-pi*cos(pi.*x));
% r = @(x)(0.*x + 1);
% g = @(x)(exp(sin(pi.*x)).*(1-pi^2*sin(pi.*x)));
% n = 64;
% 
% [f, x] = resoldreODE(n, C, p, q,r,g)
% plot(x(2: end-1), f);
% hold on;
% 
% solucio = @(x)(exp(sin(pi*x)));
% %plot(x(2: end-1), solucio(x(2: end-1)));
% %axis([0 1 0 3])
% 
% %% Practica 18
% % D,omini de x 0 a 1
% g=1;
% p = @(x)(g.*x); % En cas que sigui una funcio constant, cal que apareixi la x en el function handler
% q =@(x)(0.*x + g);
% r = @(x)(0.*x);
% n = 64;
% a=0; b=1;
% C=[0, 1, 0; 0, 1, 0];
% 
% [f, x] = resoldreODEeig(n, C, p, q,r, a, b);

%%
clear all;
clc;
close all;;
g=1;
p = @(x)(0.*x+1); % En cas que sigui una funcio constant, cal que apareixi la x en el function handler
q =@(x)(0.*x);
r = @(x)(0.*x);
n = 26;
a=0; b=1;
C=[1, 0, 0; 1, 1, 0]; 

[F, x, lamb] = resoldreODEeig(n, C, p, q,r, a, b);
disp(lamb);


plot(x(2:end-1),F(1,:))