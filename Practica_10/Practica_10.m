% x = (v,r)
% fun: (v,r) --> (dr/dt, dv/dt)
format long g;
close all;
clc;
clear all;

h = 1e-2;
time = 2;
steps = time / h;
initial = [0, 0, 1, 1]';

%Amb RK4 pur
solution = RK4(initial, h, @gravFunctionV, steps);
figure;
plot(solution(1, :), solution(2, :));

%Amb AB4
solution1 = AB4(solution(:, 1:4), h, @gravFunctionV, steps);
figure;
plot(solution1(1, :), solution1(2, :));

% Calculem el que considerem la solucio "exacte"
hext = 1e-4;
points = time / hext +1;
exactSolution = RK4(initial, hext, @gravFunctionV, points);
%exactSolution = AB4(solution(:, 1:4), hext, @gravFunctionV, steps);
rext = exactSolution(1:2, end);

rr = [];
errorsRK = [];
errorsAB = [];
ns = -4:0.01:-1;
haches = 10.^ns;
hs = [];
haches = 1e-3:0.0001:0.1;

for h = haches
    points = time / h + 1; %RK i AB ens donen el punt n+1, si volem el punt a T=2 s'ha de ficar el +1

    %{
    solution1 = RK4(initial, h, @gravFunctionV, points);
    r = solution1(:, end);
    rr = [rr r];
    errors = [errors norm(r - rext)];
    %}

    % TODO: millorar
    if floor(points) == points% Comprovar si es purament enter
        hs = [hs h];
        solution1 = RK4(initial, h, @gravFunctionV, points);
        r = solution1(1:2, end);
        errorsRK = [errorsRK norm(r - rext)];

        % AB4 començant desde les primeres 4 solucions exactes.
        solution2 = AB4(solution1(:, 1:4), h, @gravFunctionV, points);
        r2 = solution2(1:2, end);
        errorsAB = [errorsAB norm(r2 - rext)];

    end

end

figure;
loglog(hs, errorsRK)
hold on
loglog(hs, errorsAB)
hold off

% Com en el grafic logaritmic, amb el ab4 aconseguim precisio de 10-4 amb h = 0.008 i abs amb 0.025

% RK4 sera 4 vegades mes costos que el AB4

%% Section b)
% Li fixem la h = 0.0001


%{
 solutionAngle = newton2(pi/6, 1e-2, 100, @funForNewton);

h = 0.0001;
points = 2/h;
initial = [ 0; 0; sqrt(2).*cos(solutionAngle); sqrt(2).*sin(solutionAngle)];
solution1 = RK4(initial, h, @gravFunctionV, points);
figure;
plot(solution1(1, :), solution1(2, :)); 
%}
