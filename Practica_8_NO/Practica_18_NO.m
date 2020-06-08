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

%%%
% clear all;
% clc;
% close all;;
% g=1;
% p = @(x)(0.*x+1); % En cas que sigui una funcio constant, cal que apareixi la x en el function handler
% q =@(x)(0.*x);
% r = @(x)(0.*x);
% n = 26;
% a=0; b=1;
% C=[1, 0, 0; 1, 1, 0]; %fiquem les C, p, q, r, al domini de (a,b)
% 
% [F, x, lamb] = resoldreODEeig(n, C, p, q,r, a, b);
% disp(lamb);
% f=F(3,:)/abs(max(F(3,:)))
% 
% plot(x(2:end-1),f)

%%

clear all;
clc;
close all;
g=1;
p = @(x)(x); % En cas que sigui una funcio constant, cal que apareixi la x en el function handler
q =@(x)(0.*x + 1);
r = @(x)(0.*x);
n = 26;
a=0; b=1;
C=[1, 0, 1; 0, 1, 0]; %fiquem les C, p, q, r, al domini de (a,b)


[F, x, lamb] = resoldreODEeig(n, C, p, q,r, a, b);


f=F(1:3,:);
figure;
plot(x(2:end-1),f)

z = 2.*lamb(2).*sqrt(x(2:end-1))

figure;
plot(x(2:end-1),besselj(0, z))



