close all;
clear all;
format long g;

n = 30;
R = 3;
VOLTATGE = 5;


N = 2 * n - 1;
A = zeros(N);
%Ara ho fem com la practica 1:
for i = 1:1:N
    %Si es parell:
    if (mod(i, 2) == 0)
        A(i, i - 1) = 1;
        A(i, i) = -1;
        A(i, i + 1) = -1;
        %La primera i la ultima definides posteriorment
    elseif (i ~= 1 && i ~= N)
        A(i, i - 1) = -R;
        A(i, i) = 2 * R;
        A(i, i + 1) = R;
    end

end

A(1, 1) = 2 * R;
A(1, 2) = R;

v = [0 * (1:N - 2), -R, 3 * R];
A(N, :) = v;
% Muntem els termes independents

n = 30;
V = zeros(2 * n - 1, 1);
V(1) = VOLTATGE;

% Aquestes dos linies son les que solucionen!!
[x4, k] = gmres(@Afun, V, 1e-10, 100); % Se li ha de passar un funciton handler.

[P, L, U] = PLU(A);
x1 = pluSolve(L, U, P, V);

disp(x4);
disp(x1);
%disp(abs(x4 - x1));


