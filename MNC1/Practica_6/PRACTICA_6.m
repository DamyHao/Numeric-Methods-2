%% a)
clear all;
format long g;


x = [0+10^-12:1/2000:1-10^-12];
%x = [-10:1/2000:10];
plot(x, fun(x), 'b');

%% b)
% Vector z on avaluem els punts
z = [0+10^-12:1/2000:1-10^-12];
for n = 100:100:400
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

%% c)
% Per trobar el nombre de nodes necessari utilitzarem un bucle while amb un
% fucnionament similar a la biseccio per més precisió y velocitat.
clear all;
z = [0+10^-12:1/2000:1-10^-12];

salt = 1000;
n = salt; % Resolucio de la interpolacio.
prevN = 0;
maxN = 100000;
direccio = 1; % Saber si estabem augmentant o disminuint resolucio.
maxError = 10^-6;
errorInter = 1;

while n <= maxN && prevN ~= n
    n
    prevN = n;
    j=[0:1:n];
    zj = cos((j*pi)/n);
    xj = (1/2)*(1+zj);
    b = baricentrica2(z, xj, fun(xj));
    errorInter = max(abs(fun(z)- b'))
    %En cas que l'error sigui massa gran, cal augmentar resolució de les n:
    if errorInter > maxError
        if direccio == 1
            n = n + salt;
        elseif direccio == -1
            direccio = 1;
            % Divisio entera per 2:
            salt = fix(salt/2);
            n = n + salt;
        end
    % Si l'error ja es mes petit que el maxim d'error, vol dir que ens hem
    % passat amb la resolució i podem ajustar més.
    elseif errorInter < maxError
        if direccio == 1
            direccio = -1;
            % Divisio entera per 2:
            salt = fix(salt/2);
            n = n - salt;
        elseif direccio == -1
            n = n - salt;
        end
    end
end
% Deixem escrits els errors maxims i les n per comprendre el procés fet.
% Finalment hem de decidir entre n i n+1:
if errorInter <= maxError
    n
else
    n+1
end
    

        
    
