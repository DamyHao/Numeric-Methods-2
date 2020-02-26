close all;
clear all;
format long g;

A = [1, 2, 3; 2, 3, 1; 2, 5, 3; 1,2,3];
UR = qrFact(A);
qrSolve(UR);
%[Q,R] = qr(A);

f = @(x)(exp(-20 * (x + 1/4).^2) + (sin(30 * x) .* exp(-20 * (x - 1/4).^2)) / 4);

punts = [1:1:100];
disp(size(punts));


x = [-1:0.01:1];

%x = linspace(-1,1, n + 1);


casosPractica = [[14,7], [28,2]];


n = 14;
m = 7;
jota = [0:1:n+1];

xj = [-1 + 2 * jota / (n+1)];

V = fliplr(vander(xj));


% Retallem la matriu:
V = V(:,1:m+1);

efes = f(xj);



