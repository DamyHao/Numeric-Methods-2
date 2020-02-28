close all;
clear all;
format long g;

%% PRACTICA 2

%FUNCIO A ESTUDIAR:
f = @(x)(exp(-20 * (x + 1/4).^2) + (sin(30 * x) .* exp(-20 * (x - 1/4).^2)) / 4);
x = [-1:0.01:1];

%visualitzem la funcio
figure (1)
plot (x,f(x))
hold on


%% B) executem el qrSolver per trobar els coeficients de l'aproximacio de f mitjançant la matriu de Vandermonde

casos=[[14,7]; [28,14]; [28,20]; [64,30]];
ii=1;

for cas = casos'
    n = cas(1);
    m = cas(2);
    
    % Creem lel sistema
    
    jota = [0:1:n];
    xj = [-1 + 2 * jota / (n)]; %nodes equispaiats
    
    V = fliplr(vander(xj)); %Creem la Vandermonde
    V = V(:,1:m); %retallem la matriu per tal de tenir els polinomis fins a grau m
    
    efes = f(xj); %Vector de la funcio analitzada ens els nodes
    
    % El sistema sera V * x = efes, el resolem:
    
    coef = qrSolve(qrFact(V), efes');
    
    p=fliplr(vander(x)); % fem la Vandermonde pero en un seguit de punts més ampli per tal de millorar la presició 
    efesinterpolador=p(:,1:m)*coef'; % és  l'equivalent al vector de efes pero extret de la interpolacio
    
    %FUNCIO A ESTUDIAR:
    f = @(x)(exp(-20 * (x + 1/4).^2) + (sin(30 * x) .* exp(-20 * (x - 1/4).^2)) / 4);
    x = [-1:0.01:1];
    
    ii= 1 + ii;
    figure(ii)
    plot(x,f(x))
    hold on
    plot(x,efesinterpolador)
    
    
    
end


%% C) Repeat B using the Chebycheb nodes

casos=[[14,7]; [28,14]; [28,20]; [64,30]];
ii=5;

for cas = casos'
    n = cas(1);
    m = cas(2);
    
    % Creem lel sistema
    
    jota = [0:1:n];
    xj = [cos(pi * jota / (n))]; %nodes cheby
    
    V = fliplr(vander(xj)); %Creem la Vandermonde
    V = V(:,1:m); %retallem la matriu per tal de tenir els polinomis fins a grau m
    
    efes = f(xj); %Vector de la funcio analitzada ens els nodes
    
    % El sistema sera V * x = efes, el resolem:
    
    coef = qrSolve(qrFact(V), efes');
    
    p=fliplr(vander(x));
    efesinterpolador=p(:,1:m)*coef';
    
    %FUNCIO A ESTUDIAR:
    f = @(x)(exp(-20 * (x + 1/4).^2) + (sin(30 * x) .* exp(-20 * (x - 1/4).^2)) / 4);
    x = [-1:0.01:1];
    
    
    ii= 1 + ii;
    figure(ii)
    plot (x,f(x))
    hold on
    plot(x,efesinterpolador)
    
    
    
end
    

%COMENTARI FINAL: Observem que amb els nodes equispaiats en els extrems (-1
%i 1) la interpolació falla i amb els de chebycheb no, tal com vam estudiar
%que passava l'any passat.
    
    
    
    
    
    
