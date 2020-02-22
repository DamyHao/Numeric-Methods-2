close all;
clear all;
format long g;

N = 11;
R = 1;
VOLTATGE = 5;


A = zeros(N);
for i=1:1:N
    %Si es parell:
    if(mod(i,2) == 0)
        A(i,i-1) = 1;
        A(i, i) = -1;
        A(i, i+1) = -1;
        %La primera i la ultima definides posteriorment
    elseif(i ~= 1 && i ~= N)
        A(i, i-1) = -R;
        A(i, i) = 2*R;
        A(i, i+1) = R;
    end
end

A(1,1) = 2*R;
A(1,2) = R;

v=[0*(1:N-2),-R,3*R];
A(N, :) = v;

V = zeros(N,1);
V(1) = VOLTATGE;


[P, L, U] = PLU(A);
x = plusolve(L,U,P,V)


A\V
