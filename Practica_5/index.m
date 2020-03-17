%Provant...
%disp(genpath([fileparts(pwd), "\Practica_1", 'PLU.m']))
%addpath(genpath(fileparts(pwd), "\Practica_2"))
%load(fullfile('..', 'Practica_1', 'PLU.m'))
addpath("../Practica_1")

A = eye(3);
PLU(A);

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

