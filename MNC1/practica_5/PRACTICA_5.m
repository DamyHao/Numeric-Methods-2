%prova per fer un for dels tres
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
end

%% c)
clear all;
for n=[8,16,24,32]
    z = [-1:1/300:1];
    l = lebesgue(n);
    %{
    Vull posar els quatre gràfics en una sola figure, i després zones
    ampliades d'un altre gràfic amb les 4 plots combinades per observar
    millor el comportament de la funció.
for plotID=1:4
        subplot(2,2,plotID);
        plot(z,l);
        %hold on;
    end
    %}
    figure;
    plot(z,l);
    titol_plots_c = sprintf('Lebesgue per a n = %d', n);
    title(titol_plots_c);
    ylabel('Valor de la funció');
    xlabel('Punts avaluats (nodes i z_k)');
    hold on;
end
hold off;
%¿comprobar su crecimiento? Amb veure-ho a la gràfica ja n'hi ha prou?

%% d)
clear all;
for n=4:2:60
    i_exp = interpolacio_exp(n);
end
