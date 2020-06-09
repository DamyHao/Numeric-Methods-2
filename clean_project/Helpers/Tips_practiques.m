%% HELPER:

% Si es parell: 
if (mod(i, 2) == 0)

    
    
% Vandermonde: pràctica 2
% V*coef'=f --> la fem servir per trobar els coeficients del
% polinomi interpolador que fara que en els nodes el polinomi sigui exacte
% a la funcio que volem interpolar es a dir p(xj)=f(xj)
jota = 0:1:n;
xj = -1 + 2 * jota / (n); %nodes equispaiats
V = fliplr(vander(xj)); %Creem la Vandermonde
V = V(:,1:m); %retallem la matriu per tal de tenir els polinomis fins a grau m
%COMENTARI FINAL: Observem que amb els nodes equispaiats en els extrems (-1
%i 1) la interpolaciï¿½ falla i amb els de chebycheb no, tal com vam estudiar
%que passava l'any passat.



% Grafica arco-iris plus text --> practica 5:
colors = jet(length(Ts)); % Rainbow colors for the plot 
colors = fliplr(colors); %We fliplr in order to have the red colors for highest temperatures and blue for the lowest
for T = Ts
plot(xarxa, vanDerWaals(xarxa), 'Color', colors(ii, :)); %here we represent an isotherm line
text(1.1, vanDerWaals(1.1), ['T = ', num2str(T)], 'FontSize', 6); %It shows which isotherm is showing the plot


% Practica 5--> a l'apartat B fem newton de la integral d'una funcio utilitzan clenshaw
% curtis de l'any passat.



% Parctica 6: Newton exploration:
dom1 = [0, pi / 2];
dom2 = [-pi / 2, pi / 2];
factor = [dom1(2); (dom2(1) - dom2(2))];
a = [dom1(1); dom2(1)];
aleatoryTimes = 1:20;
sol = [];
alphasol = [];
for alpha = alphas
    f = @(phi)([(tan(phi(1)) - alpha * (2 * sin(phi(1)) + sin(phi(2)))), (tan(phi(2)) - 2 * alpha * (sin(phi(1)) + sin(phi(2))))]);
    % f is a function which depends on phi1 and phi2 and alpha is a
    % parameter
    for i = aleatoryTimes
        aleatory = rand(2, 1); %2 random numbers bewteen 0 and 1
        phi0 = aleatory .* factor - a; %random intital condition within the domain
        %Formula per canviar de escala i moure:
        %x*(a-b) - a
        [XK, resd, it] = newtonn(phi0, 1e-6, 100, f);
        % Comprobar que estigui dintre el domini
        if XK(1, end) > dom1(1) && XK(1, end) < dom1(2) && XK(2, end) > dom2(1) && XK(2, end) < dom2(2)
            sol = [sol, XK(:, end)];
            alphasol = [alphasol, alpha];
        end
    end
end
%Continuation step llançada a diferents branques que em torbat amb la
%newton exploration:
tol = 1e-6;
itmax = 100;
Y = [];

for it = 1:2:length(MP)% We launch the countinuation step at the 3 branches of solutions found
    y0 = MP(:, it);
    y1 = MP(:, it + 1);
    y = y1;

    while y(3) < 2 && y(3) > 0 && s > 0
        [y, iconv] = continuationStep(funAlpha, y0, y1, s, tol, itmax);

       if iconv == 1 % No hem aconseguit soluciÃ³ i ajustem s
            s = s - 0.1; % Si la s arriba a 0 desistirem i no buscarem mes solucions
        else
            y0 = y1;
            y1 = y;
            % Nomes les guardem si estan dintre el domini
            if y(1, end) >= dom1(1) && y(1, end) <= dom1(2) && y(2, end) >= dom2(1) && y(2, end) <= dom2(2)
                Y = [Y, y]; %solucions
            end
        end
        plot(Y(end, :), Y(1:2, :), 'o');
        hold on
    end
end



%Practica 8: fem plot dels 3 primers eigenvectors:
%Once we have the eigenvectors we normalize, make them positive and plot the first three of
%them
f = [];
for i=1:3
    if abs(max(F(i,:)))< abs(min(F(i,:)))
        f(i,:)=F(i,:)/min(F(i,:));
    else
        f(i,:)=F(i,:)/max(F(i,:));
    end
end
%AND THE NODES:
% Finally we find the nodes. As we have seen that using evaluatinf the 
% bessel function on  2*lamb  the result is 0, we find the nodes (x) for each 
% lambda with: lambdai = lambdaj*sqrt(x)
for i = 1:3
   for j = 1:i
       disp(strcat('Lambda ' , int2str(i) , ' nodes:'))
       nodeAt = (lamb(j)./lamb(i)).^2;
       disp(nodeAt);
   end
end



%Continuationstep by Dami en el ex2016:
s = 1;
maxIt = 1000;
iterations = 0;
Y = [];
while s > 0 && iterations <= maxIt
    [y, iconv] = continuationStep(@fapLambda, y0, y1, s, 1e-6, 100);

    if iconv == 1% No hem aconseguit soluciÃ³ i ajustem s
        s = s - 0.1; % Si la s arriba a 0 desistirem i no buscarem mes solucions
    else
        y0 = y1;
        y1 = y;
        Y = [Y, y]; %solucions
    end
    plot(y(end), max(y(1:end-1)), 'o');
    hold on;
    iterations = iterations + 1;
end