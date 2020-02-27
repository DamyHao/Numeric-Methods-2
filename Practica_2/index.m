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
    V = V(:,1:m+1); %retallem la matriu per tal de tenir els polinomis fins a grau m
    
    efes = f(xj); %Vector de la funcio analitzada ens els nodes
    
    % El sistema sera V * x = efes, el resolem:
    
    x= qrSolve(qrFact(V), efes');
    
    y=polyval(x,xj); %polinomi interpolador
    
    
    ii= 1 + ii;
    figure(ii)
    plot(xj,y)
   
    
    
end

    
    
    
    
    
    
    
    
    
