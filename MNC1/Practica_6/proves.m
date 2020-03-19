clear all;
close all
format long g;
z = [0+10^-12:1/2000:1-10^-12];
%{
for n = 100:100:400
    j=[0:1:n];
    zj = cos((j*pi)/n);
    xj = (1/2)*(1+zj);
    P = MLagrange(z, xj);
    fun_transposada = fun(xj)';
    fun_interpolada = P * fun_transposada;
    figure;
    plot(z, fun_interpolada);
    titol_plot_int = sprintf('FunciÃ³ interpolÂ·lada amb Lagrange per a n = %d', n);
    title(titol_plot_int);
    xlabel('z');
    ylabel('funciÃ³');
    hold on;
    error = abs(fun(z)- fun_interpolada');
    figure;
    semilogy(z,error,'k');
    titol_plot_err = sprintf('Error local amb Lagrange per a n = %d', n);
    title(titol_plot_err);
    xlabel('z');
    ylabel('E_k')
    hold on;
end
hold off;

%}


z = [0+10^-12:1/2000:1-10^-12];
for n = 403:100:403
    j=[0:1:n];
    zj = cos((j*pi)/n);
    xj = (1/2)*(1+zj);
    b = baricentrica2(z, xj, fun(xj));
    figure;
    plot(z, b);
    titol_plot_int = sprintf('Funció interpola·lada amb Lagrange per a n = %d', n);
    title(titol_plot_int);
    xlabel('z');
    ylabel('funció');
    hold on;
    error = abs(fun(z)- b');
    figure;
    semilogy(z,error,'k');
    titol_plot_err = sprintf('Error local amb Lagrange per a n = %d', n);
    title(titol_plot_err);
    xlabel('z');
    ylabel('E_k')
    hold on;
end
hold off;
maxError= max(error);