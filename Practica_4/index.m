close all;
clear all;
format long g;

n = 30;
R = 3;
VOLTATGE = 5;

addpath('../Practica_1'); % Podem utilitzar totes les funcions de la practica 1.

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

V = zeros(N, 1);
V(1) = VOLTATGE;
%TODO: EL SISTEMA NO ES RESOL SI LI POSEM RESISTENCIA 3
[P, L, U] = PLU(A);

x = pluSolve(L, U, P, V);

%Comprovem que el Afun fa el mateix que la matriu A
vectorDeProva = ones(N,1);

disp('Comprobem que fan el mateix Afun i Ax')
disp(norm(Afun(vectorDeProva)-A*vectorDeProva))

% Aquestes dos linies son les que solucionen!!
[x4, k] = gmres(@Afun, V, 1e-15, 100); % Se li ha de passar un funciton handler.
[x5, k] = mygmres(@Afun, V, 1e-15, 100); % Se li ha de passar un funciton handler.

% Mirem les diferencies entre el resultat de la practica 1 i la 4
disp(norm(x4 - x'));
disp(norm(x'- A\V));
disp(norm(x4- A\V));
disp(norm(x5- A\V));

solucioReal = A\V;
y = -2:-1:-10;
tolerancies = 10.^y;

ks = [];
ks2 = [];
for tol = tolerancies
    [x4, ktemp] =  gmres(@Afun, V, tol, 100);
    [x5, ktemp2] = mygmres(@Afun, V, tol, 100); % Se li ha de passar un funciton handler.
    
    ks = [ks, ktemp];
    ks2 = [ks2, ktemp2];
end

figure();
semilogx(tolerancies, ks);
ylabel('Dimensio de krylov');
xlabel("Tolerancia");

figure();
semilogx(tolerancies, ks2);
ylabel('Dimensio de krylov');
xlabel("Tolerancia");