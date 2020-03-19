%% a)
close all;
clear all;
n = 10;
ncont= [0:1:n];
x = -1+(2*ncont)/n;
z = [-1:1/300:1];
P = MLagrange(z, x);

%% b)
clear all;
for n=[3,6,9]
    ncont= [0:1:n];
    x = -1+(2*ncont)/n;
    z = [-1:1/300:1];
    P = MLagrange(z, x);
    
    figure;
    for columna = 1:1:(n+1)
        plot(z, P(:,columna));
        hold on;
    end
    titol_plots_b = sprintf('Lagrange per a n = %d', n);
    title(titol_plots_b);
    ylabel('Polinomis cardinals de Lagrange');
    xlabel('punts avaluats (nodes i zk)');
    hold off;
end;

%% c)
clear all
figure;
for n = [8, 16, 24, 32]
    ncont= [0:1:n];
    x = -1+(2*ncont)/n;
    z = [-1:1/300:1];
    P = MLagrange(z, x);
    y = [];
    for i=1:1:length(z)
        fila = abs(P(i,:));
        y = [y sum(fila)];
    end;
    plot(z, y);
    hold on;
end
hold off;

%% d)
clear all
% Ja tenim feta la matriu dels polinomis cardinals de Lagrange, ara falta
% fer el vector amb els valors de l'exponencial.
figure('Name','Exponencial');
maxError = [];
% A fora per optimitzar
z = [-1:1/300:1];
expZ = exp(z)';

for n = 4:2:60
    ncont= [0:1:n];
    x = -1+(2*ncont)/n;
    P = MLagrange(z, x);
    % Vector exponencial en x:
    expX = exp(x);
    expX = expX';
    P = P * expX;
    y = [];
    for i=1:1:length(z)
        % Com el exercici anterior pero treiem el valor absolut.
        % El vector y resulta ser les avaluacions del polinomi
        % interpolador.
        fila = P(i,:);
        y = [y sum(fila)];
    end;
    if n==4
        plot (z, y, 'r');
    else
        plot(z, y);
    end;
    hold on;
    
    % Aqui calculem el error maxim:
    maxError = [maxError max(abs(expZ' - y))];
end
hold off;

figure;
n = [4:2:60];
plot(n, maxError);
%hold on;
figure;
n = n';
twos = 2.*ones(length(n),1);
probarLey = (twos.^n)/(n.*log(n));
% Perquè probarLey es un matriu 29x29????
plot(n, probarLey(:,length(n)), 'r');



