%% Damia Casas Casajuana P2
clear all;
close all;
format long g;

%% a)
xi = [1:1/100:20];
   
for alfa = [3 7 10]  
    figure();
 
    yleft = xi.^2 + 4.*xi.*sin(xi) + (2*sin(xi)).^2;
    plot(xi, yleft);
    hold on;
    
    yright = alfa*xi;
    plot(xi, yright);
    hold off;
end

%% b)
% Si hem de trobar 8 xifres signfigicatives, augmento la resolució dels
% punt fins el maxim que em permet el Matlab:

xmax = 20;
xmin = 1;
xi = [xmin:1/10000000:xmax];
itmax = 100;
tol = 10^-8;
% Matriu de solucions en que es guardaran les diferents solucions de totes
% les alfes (les alfes seran les files i les solucoins estaran en columnes
MSolucions = [];
solucio = [];
h= 10.^-10;

for alfa = [3 7 10]
    solucions = [];
    
    % Defineixo la funcio respecte x:
    funcio = @(x)((x.^2 + 4.*x.*sin(x) + (2*sin(x)).^2)-(alfa*x));
    
    % M'asseguro de llançar Newton en un lloc on hi ha solució segons
    % Bolzano:
    xa = 1;
    xb = xa + 0.25;
    
    while (xb < xmax + 0.26)
        xa = xb;
        xb = xb + 0.25;
        % Si compleix els requisits:
        if funcio(xa)*funcio(xb) < 0
            % Llanço Newton:
            it=0;
            tolk=1;
            xk = [(xa+xb)/2];
           
            while (it<itmax && tolk>tol)
                % Primer calculo la derivada en el punt correcte
                f1 = funcio(xk(end)+h);
                f2 = (funcio(xk(end)));
                d = (f1-f2)/h;
                xk = [xk xk(end)-(funcio(xk(end))/d)];
                tolk = abs(xk(end-1)-xk(end));
                it = it+1;
            end
            % D'aquesta manera automaticament surten les solucions en el
            % publish i per la alfa que s'ha calculat.
           alfa
           xk(end)
           % solucio = [solucio + xk(end)];
        end 
    end
    % MSolucions(alfa) = solucio';
end

%% c)

