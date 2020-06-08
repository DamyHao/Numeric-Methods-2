%%a)
clear all;
close all
[D, x] = chebdiff(16, -1, 1);

a = (1-x.^2);
A = diag(a);

b = (-4.*x);
B = diag(b);

c = 2;
C = -2*eye(17);
%C = diag(-2)

d = 1-x.^2;
DD = diag(d);

global L;
L = (A*(D^2) + B*D + C);


% B No es invertible per tant:

LL = inv(DD(2:end-1, 2:end-1))*L(2:end-1, 2:end-1);
[EVEC, EVAL] = eig(LL);

lamb = diag(EVAL); % Lamb es un vector perq diag funciona en les dos direccions de conversio.
[foo,ii] = sort(lamb) ; lamb = lamb(ii) ; EVEC = EVEC(:,ii);



%% b
x0 = [1;1];
[XK, conv, it] = newtonnConv(x0, 1e-8, 15, @funForNewton);




%% c)
close all;
clc
p = @(x)(0 .* x + 1); % En cas que sigui una funcio constant, cal que apareixi la x en el function handler
q = @(x)(0 .* x + 1);
r = @(x)(0 .* x +0);
g = @(x)(0 .*x + 0);
C= [1, 0, 0; 1, 0, 0];
n = 24;
[f, x] = resoldreODE(n, C, p, q, r, g, -1, 1);

plot(x(2:end-1), f);
n = 126;
initial = ones(n+1, 1).*0.5;
global D nodes
[D, x] = chebdiff(n, -1, 1);
D = D^2;
%[XK2, conv2, it2] = newtonnConv(initial, 1e-14, 500, @fun2);
nodes = exp(x);
%plot(x, XK2(:, end));

%XK2(14, end)

[XK3, conv3, it3] = newtonnConv(initial, 1e-14, 500, @fun3);
plot(x, XK3(:, end));

plot(x, x+XK3(:, end));
AA = x+XK3(:, end);
function out = fun2 (fs)
global D
out = [fs(1); D*fs-exp(fs); fs(end)];
end

function out = fun3 (gs)
global D nodes
out = [gs(1); D*gs-exp(gs).*nodes; gs(end)];
end

function result = funForNewton(position)

force = [1.2; -0.9];
% x ; y
restForce = [1; -1];
actualForce = restForce + (1/norm(position)^3)*position;
result = force-actualForce;
end
