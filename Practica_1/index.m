close all;
clear all;
format long g;

N = 11;
R = 1;
VOLTATGE = 5;

%% a)

for n = 5:5:30
    N = 2 * n - 1;
    A = zeros(N);
    % Muntem la matriu (les matrius que muntem han de ser amb files imparells sempre)
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

    % Aquestes dos linies son les que solucionen!!
    [P, L, U] = PLU(A);
    x = pluSolve(L, U, P, V);

    if (n == 10)
        figure()
        plot(1:1:N, x)%Si que s'observa una disminucio de la intensitat
        ylabel('Intensitat')
        xlabel('k')
    end

end

%% b)
% Provem amb una N molt gran:
N = 1001; % Imparell
A = zeros(N);
% Muntem la matriu
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

% Aquestes dos linies son les que solucionen!!
[P, L, U] = PLU(A);
x = pluSolve(L, U, P, V);

figure()
loglog(1:1:N, x)%Aquet plot no caldria, ya que nomes fem servir I(1)
ylabel('Intensitat')
xlabel('k')

% Veiem que la primera intensitat es 1.
disp(x(1));
disp(VOLTATGE/x(1));